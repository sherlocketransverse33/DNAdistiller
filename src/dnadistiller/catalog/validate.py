"""Standalone catalog validation, run in CI and usable by contributors.

    uv run python -m dnadistiller.catalog.validate

`load_catalog` already enforces the schema, so this adds the checks that are too
opinionated to block a library import but should still block a merge: duplicate
rsIDs across markers, interpretations that do not cover the plausible genotypes,
and citations that are not in a resolvable form.

A typo'd rsID is the reason this exists. It does not raise, it does not warn, and
it does not show up in any test that uses a fixture. It just reports "not tested"
to every user forever.
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict

from ..models import Evidence, Marker
from . import CatalogError, load_catalog

#: PMID:12345678 or a DOI. Anything else is probably a URL that will rot, or a
#: paper title that cannot be checked mechanically.
_CITATION_RE = re.compile(r"^(PMID:\d{1,9}|DOI:10\.\d{4,9}/\S+|10\.\d{4,9}/\S+)$", re.IGNORECASE)

_VALID_ALLELES = set("ACGTID-0")


def check_citations(marker: Marker) -> list[str]:
    problems = []
    for citation in marker.citations:
        if not _CITATION_RE.match(citation.strip()):
            problems.append(
                f"citation {citation!r} is not a PMID or DOI. "
                "Use 'PMID:12345678' or '10.1000/xyz' so it can be resolved."
            )
    return problems


def check_interpretation_coverage(marker: Marker) -> list[str]:
    """Warn when a marker does not cover all three genotypes at a bi-allelic site.

    A catalog entry that only describes the risk genotype produces a profile in
    which every line sounds like bad news, because the common result silently
    renders as UNKNOWN_GENOTYPE instead of "typical".
    """
    if marker.derivation:
        return []  # derived markers define their own key space

    alleles = {a for genotype in marker.interpretations for a in genotype if a in "ACGT"}
    if len(alleles) != 2:
        return []

    first, second = sorted(alleles)
    expected = {first + first, first + second, second + second}
    missing = expected - set(marker.interpretations)
    if missing:
        return [
            f"interpretations cover {sorted(marker.interpretations)} but not "
            f"{sorted(missing)}. A user with an uncovered genotype gets no result. "
            "Describe the common genotype too, even if it is unremarkable."
        ]
    return []


def check_effect_allele(marker: Marker) -> list[str]:
    if marker.effect_allele is None:
        return []
    present = {a for genotype in marker.interpretations for a in genotype}
    if marker.effect_allele not in present:
        return [
            f"effect_allele {marker.effect_allele!r} does not appear in any "
            f"interpretation key. Either it is on the wrong strand, or it is a typo. "
            "Copy counts will all read as zero."
        ]
    return []


def check_marker_is_not_inverted(marker: Marker) -> list[str]:
    """Catch a marker whose alleles have been transcribed the wrong way round.

    This check exists because of a real, shipped, published failure in another
    consumer-DNA tool. It labelled TPMT rs1142345 `TT` as "poor metaboliser,
    standard thiopurine doses can cause fatal myelosuppression". `T` is the
    reference allele at that position, carried by roughly 96% of people. The
    interpretation was inverted, so almost every user who ran it was told they
    were at risk of a fatal drug reaction. It reached their sample report and
    their test suite asserted the wrong tier, which is what let it survive.

    The reason a strand fallback cannot save you here is worth stating plainly,
    because it is counterintuitive: complement matching makes the catalog's
    orientation stop mattering *for matching* while it still entirely determines
    *the meaning attached*. A minus-strand file complements neatly onto the
    inverted entry and gets the same wrong answer. Strand tolerance hides
    orientation errors rather than surfacing them.

    So the check is on frequency, not on strand. An effect allele whose
    homozygote is carried by most of the population is almost always a
    transcription inversion: real effect alleles are usually the minor allele,
    and a "risk variant" that 96% of people carry is a contradiction in terms.
    Where a marker genuinely has a common effect allele, state the frequency and
    the reviewer will see this fire and can say so in the pull request.
    """
    if marker.effect_allele is None:
        return []

    homozygous_effect = marker.effect_allele * 2
    interpretation = marker.interpretations.get(homozygous_effect)
    if interpretation is None or interpretation.population_frequency is None:
        return []

    if interpretation.population_frequency > 0.5:
        return [
            f"genotype {homozygous_effect} carries the effect allele twice and is "
            f"declared at {interpretation.population_frequency:.0%} population "
            "frequency. An effect allele carried by most people usually means the "
            "alleles are inverted: check which allele is the reference at this "
            "position on dbSNP. If the effect allele really is the common one, say "
            "so in the pull request so a reviewer can confirm it."
        ]
    return []


def check_frequencies_are_plausible(marker: Marker) -> list[str]:
    """Population frequencies across a marker's genotypes should roughly sum to 1.

    Loose on purpose: the numbers come from different papers and different
    populations, so they will not add up exactly. What this catches is the
    order-of-magnitude slip, such as a percentage written where a proportion was
    meant, which is also the shape of an inversion.
    """
    frequencies = [
        i.population_frequency
        for i in marker.interpretations.values()
        if i.population_frequency is not None
    ]
    if len(frequencies) < 2:
        return []

    total = sum(frequencies)
    if not 0.7 <= total <= 1.3:
        return [
            f"population frequencies sum to {total:.2f} across {len(frequencies)} "
            "genotypes. They should roughly total 1. Check for a percentage written "
            "as a proportion, or a missing genotype."
        ]
    return []


def check_alleles(marker: Marker) -> list[str]:
    if marker.derivation:
        return []
    problems = []
    for genotype in marker.interpretations:
        stray = set(genotype.upper()) - _VALID_ALLELES
        if stray:
            problems.append(f"interpretation key {genotype!r} contains {sorted(stray)}")
    return problems


def check_strong_claims(marker: Marker) -> list[str]:
    """Hold the strong grade to a higher bar than the loader does."""
    problems = []
    if marker.evidence is Evidence.STRONG and not marker.effect_size:
        problems.append(
            "graded strong but states no effect_size. A strong claim should say how "
            "large the effect is, so a reader can judge whether it matters to them."
        )
    if marker.evidence in (Evidence.STRONG, Evidence.MODERATE) and not marker.ancestry_note:
        problems.append(
            "states no ancestry_note. Most GWAS is European-ancestry-biased and "
            "transfers poorly; say which populations this came from."
        )
    return problems


def check_duplicate_rsids(markers: list[Marker]) -> list[str]:
    """Two markers reading the same rsID is usually a copy-paste, not a design."""
    owners: dict[str, list[str]] = defaultdict(list)
    for marker in markers:
        for rsid in marker.rsids:
            owners[rsid.lower()].append(marker.id)

    return [
        f"{rsid} is claimed by more than one marker: {', '.join(ids)}. "
        "If that is deliberate, say so in a comment."
        for rsid, ids in sorted(owners.items())
        if len(ids) > 1
    ]


PER_MARKER_CHECKS = (
    check_citations,
    check_interpretation_coverage,
    check_effect_allele,
    check_marker_is_not_inverted,
    check_frequencies_are_plausible,
    check_alleles,
    check_strong_claims,
)


def main() -> int:
    try:
        catalog = load_catalog()
    except CatalogError as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 1

    problems: list[str] = []

    for marker in catalog:
        for check in PER_MARKER_CHECKS:
            problems.extend(f"{marker.id}: {problem}" for problem in check(marker))

    problems.extend(check_duplicate_rsids(catalog))

    if problems:
        print(f"FAIL: {len(problems)} problem(s) in {len(catalog)} markers\n", file=sys.stderr)
        for problem in problems:
            print(f"  - {problem}", file=sys.stderr)
        return 1

    by_evidence: dict[Evidence, int] = defaultdict(int)
    for marker in catalog:
        by_evidence[marker.evidence] += 1

    print(f"OK: {len(catalog)} markers")
    for evidence in Evidence:
        print(f"  {evidence.value:>8}: {by_evidence[evidence]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
