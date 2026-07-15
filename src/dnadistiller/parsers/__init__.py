"""Format detection and dispatch for consumer DNA exports.

Providers all ship a variant of "rsID, chromosome, position, call" and differ in
every detail around it: separators, whether the header is commented, whether
alleles arrive as one field or two, and what a failed call looks like. The
parsers normalise all of that into `Genome`.

Detection sniffs content rather than trusting the filename, because by the time
a file reaches us it has usually been renamed, and `genome.txt` says nothing.
Guessing wrong is not a cosmetic failure — misreading which column holds the
genotype produces a profile that is wrong rather than empty — so
`detect_format` refuses rather than falls back to a default.
"""

from __future__ import annotations

from pathlib import Path

from ..models import Genome
from .providers import (
    AncestryDNAParser,
    FamilyTreeDNAParser,
    GenericParser,
    LivingDNAParser,
    MyHeritageParser,
    Parser,
    TwentyThreeAndMeParser,
    VCFParser,
    diagnose,
)

#: Order matters. 23andMe and AncestryDNA both use tab-separated data under
#: comment headers, and AncestryDNA is distinguished only by its uncommented
#: column-header row and split allele columns — so it must be offered the file
#: before the more permissive 23andMe parser claims it.
PARSERS: tuple[type[Parser], ...] = (
    VCFParser,
    AncestryDNAParser,
    TwentyThreeAndMeParser,
    MyHeritageParser,
    FamilyTreeDNAParser,
    LivingDNAParser,
    # Always last. It matches on structure rather than on a vendor name, so it
    # would happily claim a file that a specific parser understands better.
    GenericParser,
)


class UnknownFormatError(Exception):
    """Raised when no parser recognises a file.

    Carries the head of the file, because the useful next step is almost always
    for a human to look at the first few lines — and for us to add a parser.
    """


def _check_not_an_archive(path: Path) -> None:
    """Catch the still-zipped download before it becomes a confusing parse error.

    Providers hand out a .zip, and people rename it or their browser does. A
    file called `genome.txt` whose first bytes are `PK` is common enough to be
    worth naming, since the alternative is an unrecognised-format error that
    sends someone hunting for a bug in their own file.
    """
    with path.open("rb") as handle:
        magic = handle.read(4)

    if magic[:2] == b"PK":
        raise UnknownFormatError(
            f"{path.name} is a ZIP archive, not a raw DNA file.\n\n"
            f"Unzip it first:  unzip {path.name}\n"
            "Then run dnadistiller against the .txt inside."
        )
    if magic[:2] == b"\x1f\x8b":
        raise UnknownFormatError(
            f"{path.name} is gzip-compressed, not a raw DNA file.\n\n"
            f"Decompress it first:  gunzip {path.name}"
        )


def sniff(path: Path, lines: int = 40) -> list[str]:
    """Read the first `lines` lines for detection, tolerating binary junk.

    `errors="replace"` rather than strict: some exports carry a BOM or a stray
    non-UTF-8 byte in the header comments, and refusing to even look at a file
    over one bad byte in a comment would be obnoxious.
    """
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        # strict=False: the file is expected to outrun the range. That is the point.
        return [line.rstrip("\n") for _, line in zip(range(lines), handle, strict=False)]


def detect_format(path: Path) -> type[Parser]:
    """Identify which parser handles this file.

    Raises UnknownFormatError rather than guessing. A wrong guess yields
    confidently wrong genetics; an error yields a bug report.
    """
    _check_not_an_archive(path)

    head = sniff(path)
    for parser in PARSERS:
        if parser.matches(head):
            return parser

    preview = "\n".join(head[:8]) or "(file is empty)"
    raise UnknownFormatError(
        f"Could not recognise the format of {path.name}.\n\n"
        f"First lines:\n{preview}\n\n"
        f"Supported: {', '.join(p.provider for p in PARSERS)}.\n"
        "If this is a raw export from a provider we do not handle yet, please open an "
        "issue. Include the first 10 lines with your rsIDs removed, not the file itself."
    )


def parse(
    path: Path,
    parser: type[Parser] | None = None,
    *,
    expect_full_export: bool = True,
) -> Genome:
    """Parse a raw DNA export into a Genome, detecting the format if not given.

    Diagnostics run here rather than inside each parser, so every provider gets
    the same scrutiny and a new parser cannot forget to be scrutinised. The
    result is attached to the Genome instead of printed, because the same facts
    have to reach a terminal, a profile, and whatever comes next.

    `expect_full_export=False` suppresses the truncated-file warning for
    fixtures and fragments, where being small is the point.
    """
    chosen = parser or detect_format(path)
    genome = chosen().parse(path)
    genome.issues = diagnose(genome, expect_full_export=expect_full_export)
    return genome


__all__ = [
    "PARSERS",
    "Parser",
    "UnknownFormatError",
    "detect_format",
    "diagnose",
    "parse",
    "sniff",
]
