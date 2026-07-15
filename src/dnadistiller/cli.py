"""Command-line interface.

Defaults here are the product. Most people will run `dnadistiller profile my.txt`
once, take whatever it prints, and paste it somewhere — so the default has to be
the choice we would defend on their behalf: STANDARD tier, no sensitive markers,
moderate-or-better evidence only. Every riskier option exists, and every one is
something you have to ask for by name.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console
from rich.table import Table

from . import __version__
from .catalog import CatalogError, load_catalog
from .models import Build, Category, Evidence, Genome, Severity, Tier
from .output import RENDERERS
from .parsers import UnknownFormatError, detect_format, parse
from .profile import build_profile, select_markers

app = typer.Typer(
    name="dnadistiller",
    help="Build a longevity-focused genetic profile from your DNA file, locally.",
    add_completion=False,
    no_args_is_help=True,
)

# stderr, so that `dnadistiller profile x.txt > profile.md` redirects the profile and
# leaves warnings on the terminal where they can still be read.
console = Console(stderr=True)
stdout = Console()

_BUILD_LABELS: dict[Build, str] = {
    Build.NCBI36: "[yellow]NCBI36 / hg18 (very old export)[/yellow]",
    Build.GRCH37: "GRCh37 / hg19",
    Build.GRCH38: "GRCh38",
    Build.UNKNOWN: "[yellow]not stated in file[/yellow]",
}


def _version_callback(value: bool) -> None:
    if value:
        stdout.print(f"dnadistiller {__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: Annotated[
        bool,
        typer.Option("--version", callback=_version_callback, is_eager=True, help="Show version."),
    ] = False,
) -> None:
    """dnadistiller reads your DNA file locally. It has no network code."""


# ---------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------


@app.command()
def check(
    path: Annotated[Path, typer.Argument(help="Your raw DNA export file.", exists=True)],
) -> None:
    """Verify a DNA file parses correctly, without interpreting anything.

    Worth running first. It answers "did this tool understand my file at all",
    which is a different question from "what does my file say" — and if the
    answer is no, every result downstream would be wrong in ways that are hard
    to notice.
    """
    try:
        parser = detect_format(path)
    except UnknownFormatError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(code=1) from None

    console.print(f"Detected format: [bold]{parser.provider}[/bold]")

    genome = parse(path, parser)

    if not len(genome):
        console.print(
            "[red]Parsed 0 variants. The file matched a format but contained no data.[/red]"
        )
        raise typer.Exit(code=1)

    table = Table(title=f"{path.name}", show_header=False, box=None)
    table.add_row("Provider", parser.provider)
    table.add_row("Variants", f"{len(genome):,}")
    table.add_row("Reference build", _BUILD_LABELS[genome.build])

    no_calls = sum(1 for g in genome if g.is_no_call)
    rate = no_calls / len(genome)
    colour = "green" if rate < 0.02 else "yellow" if rate < 0.05 else "red"
    table.add_row("No-calls", f"[{colour}]{no_calls:,} ({rate:.1%})[/{colour}]")
    stdout.print(table)

    counts = Counter(g.chromosome for g in genome)
    chrom_table = Table(title="Variants per chromosome")
    chrom_table.add_column("Chr")
    chrom_table.add_column("Count", justify="right")
    for chromosome in sorted(counts, key=_chromosome_sort_key):
        chrom_table.add_row(chromosome, f"{counts[chromosome]:,}")
    stdout.print(chrom_table)

    _report_issues(genome)
    _report_coverage(genome)
    _report_blind_spots(genome)


_SEVERITY_STYLE: dict[Severity, str] = {
    Severity.ERROR: "red",
    Severity.WARNING: "yellow",
    Severity.INFO: "dim",
}


def _report_issues(genome: Genome) -> None:
    """Print parse diagnostics worst-first.

    Ordered by severity rather than by discovery, because the one that means
    "do not trust anything below this" has to be the one you see first.
    """
    if not genome.issues:
        console.print("\n[green]No parse issues.[/green]")
        return

    console.print()
    for issue in sorted(genome.issues, key=lambda i: list(Severity).index(i.severity)):
        style = _SEVERITY_STYLE[issue.severity]
        console.print(
            f"[{style}]{issue.severity.value.upper()}[/{style}] [dim]({issue.code})[/dim] {issue.message}"
        )


def _report_blind_spots(genome: Genome) -> None:
    """State what this provider's array structurally cannot see.

    Printed on every run, including clean ones. These are not a failure report:
    they are the difference between "we found nothing" and "there is nothing",
    and the user cannot tell those apart without being told.
    """
    if not genome.blind_spots:
        return
    # "What this X file" rather than "What a X file": the provider name is data,
    # so the article would have to agree with it, and "a AncestryDNA" is what you
    # get when it does not.
    console.print(f"\n[bold]What this {genome.source} file cannot tell you[/bold]")
    for spot in genome.blind_spots:
        console.print(f"[dim]  - {spot}[/dim]")


def _chromosome_sort_key(chromosome: str) -> tuple[int, str]:
    """Sort 1-22 numerically, then X, Y, MT — rather than 1, 10, 11, 2."""
    order = {"X": 23, "Y": 24, "MT": 25}
    if chromosome.isdigit():
        return (int(chromosome), "")
    return (order.get(chromosome, 99), chromosome)


def _report_coverage(genome: Genome) -> None:
    """Say how much of the catalog this particular chip can actually answer.

    Chip content varies by provider and year, so coverage is a property of the
    user's file rather than of the tool. Reporting it up front avoids the
    reasonable-but-wrong reading of a sparse profile as a clean bill of health.
    """
    try:
        catalog = load_catalog()
    except CatalogError:
        return

    covered = sum(1 for m in catalog if any(rsid in genome for rsid in m.rsids))
    console.print(
        f"\nCatalog coverage: [bold]{covered}[/bold] of {len(catalog)} markers are present "
        f"on this chip."
    )
    if covered < len(catalog):
        console.print(
            "[dim]The rest were never tested by this provider. That is normal: arrays "
            "test preselected positions, they do not sequence.[/dim]"
        )


# ---------------------------------------------------------------------------
# profile
# ---------------------------------------------------------------------------


@app.command()
def profile(
    path: Annotated[Path, typer.Argument(help="Your raw DNA export file.", exists=True)],
    tier: Annotated[
        Tier,
        typer.Option(
            "--tier",
            "-t",
            help="How much detail to disclose. minimal=interpretations only, "
            "standard=genes and copy counts, full=rsIDs and genotypes.",
        ),
    ] = Tier.STANDARD,
    categories: Annotated[
        list[Category] | None,
        typer.Option("--category", "-c", help="Limit to these topics. Repeatable. Default: all."),
    ] = None,
    include_apoe: Annotated[
        bool,
        typer.Option(
            "--include-sensitive",
            help="Include markers with results that are hard to un-know, notably APOE "
            "and Alzheimer's risk. Off by default.",
        ),
    ] = False,
    min_evidence: Annotated[
        Evidence,
        typer.Option("--min-evidence", help="Drop markers with weaker support than this."),
    ] = Evidence.MODERATE,
    output_format: Annotated[
        str, typer.Option("--format", "-f", help="Output format: md or json.")
    ] = "md",
    out: Annotated[
        Path | None,
        typer.Option("--out", "-o", help="Write to a file instead of stdout."),
    ] = None,
    no_prompts: Annotated[
        bool, typer.Option("--no-prompts", help="Omit the suggested follow-up questions.")
    ] = False,
) -> None:
    """Build a shareable longevity profile from your DNA file."""
    renderer = RENDERERS.get(output_format.lower())
    if renderer is None:
        console.print(
            f"[red]Unknown format {output_format!r}. Choose from: {', '.join(sorted(RENDERERS))}[/red]"
        )
        raise typer.Exit(code=1)

    try:
        catalog = load_catalog()
    except CatalogError as exc:
        console.print(f"[red]Catalog error: {exc}[/red]")
        raise typer.Exit(code=1) from None

    try:
        genome = parse(path)
    except UnknownFormatError as exc:
        console.print(f"[red]{exc}[/red]")
        raise typer.Exit(code=1) from None

    selected = select_markers(
        catalog,
        categories=set(categories) if categories else None,
        min_evidence=min_evidence,
        include_sensitive=include_apoe,
    )

    if not selected:
        console.print("[yellow]No markers matched those filters. Nothing to report.[/yellow]")
        raise typer.Exit(code=1)

    result = build_profile(genome, selected, tier=tier)
    rendered = renderer(result, include_prompts=not no_prompts)

    if out:
        out.write_text(rendered, encoding="utf-8")
        console.print(f"Wrote {out} ({len(result.reportable)} findings).")
        _warn_about_sharing(tier, out)
    else:
        stdout.print(rendered, end="", highlight=False, markup=False)
        _warn_about_sharing(tier, None)


def _warn_about_sharing(tier: Tier, out: Path | None) -> None:
    """Say plainly what this file is before anyone pastes it anywhere.

    Placed after the output rather than before it, because a warning printed
    before a wall of Markdown has scrolled away by the time the decision to
    paste gets made.
    """
    console.print()
    if out:
        console.print(
            f"[bold]{out}[/bold] contains your genetic results. Review it before sharing."
        )

    if tier is Tier.FULL:
        console.print(
            "[yellow]This is a FULL profile: it contains rsIDs and genotypes. Anyone "
            "holding it can look up every variant, and enough genotypes together can "
            "single you out. Share it with a person or service you would trust with a "
            "medical record.[/yellow]"
        )
    elif tier is Tier.STANDARD:
        console.print(
            "[dim]This is a STANDARD profile: gene names and copy counts, no raw "
            "genotypes. Less identifying than FULL, but not anonymous. "
            "Try --tier minimal to share less.[/dim]"
        )
    else:
        console.print(
            "[dim]This is a MINIMAL profile: interpretations only. This is the least "
            "you can share and still have a useful conversation, though it is still "
            "information about your body, and it is still not anonymous.[/dim]"
        )


# ---------------------------------------------------------------------------
# markers
# ---------------------------------------------------------------------------


@app.command()
def markers(
    category: Annotated[
        Category | None, typer.Option("--category", "-c", help="Limit to one topic.")
    ] = None,
    show_weak: Annotated[
        bool, typer.Option("--show-weak", help="Include markers graded weak.")
    ] = False,
) -> None:
    """List the marker catalog and the evidence behind each entry.

    Reads no DNA file. This is the tool's claims laid out for inspection before
    you decide whether to trust it with anything.
    """
    try:
        catalog = load_catalog()
    except CatalogError as exc:
        console.print(f"[red]Catalog error: {exc}[/red]")
        raise typer.Exit(code=1) from None

    shown = [
        m
        for m in catalog
        if (category is None or m.category is category)
        and (show_weak or m.evidence is not Evidence.WEAK)
    ]

    table = Table(title=f"dnadistiller marker catalog ({len(shown)} of {len(catalog)})")
    table.add_column("Gene", style="bold")
    table.add_column("Trait")
    table.add_column("Category")
    table.add_column("Evidence")
    table.add_column("rsID", style="dim")

    colours = {Evidence.STRONG: "green", Evidence.MODERATE: "yellow", Evidence.WEAK: "red"}
    for marker in shown:
        flag = " [red]*[/red]" if marker.sensitive else ""
        table.add_row(
            marker.gene + flag,
            marker.name,
            str(marker.category),
            f"[{colours[marker.evidence]}]{marker.evidence}[/{colours[marker.evidence]}]",
            ", ".join(marker.rsids),
        )

    stdout.print(table)
    if any(m.sensitive for m in shown):
        console.print("[dim]* Sensitive: excluded unless --include-sensitive is passed.[/dim]")
    if not show_weak:
        weak = sum(1 for m in catalog if m.evidence is Evidence.WEAK)
        if weak:
            console.print(
                f"[dim]{weak} marker(s) graded weak are hidden. --show-weak to see them "
                "and why they are graded that way.[/dim]"
            )


if __name__ == "__main__":  # pragma: no cover
    sys.exit(app())
