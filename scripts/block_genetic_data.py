#!/usr/bin/env python3
"""Pre-commit hook: refuse to commit anything that looks like real genetic data.

`.gitignore` stops the accident where a DNA file sits in the working tree and
gets swept up by `git add .`. It does nothing about `git add -f`, a file placed
inside `tests/fixtures/` where the ignore rules are relaxed, or a path nobody
predicted. This hook looks at content instead of names.

It is deliberately noisy in one direction. A false positive costs a contributor
one `--no-verify` and thirty seconds of annoyance. A false negative puts a real
person's genome, and their relatives' genomes, into a public git history that
cannot be rewritten out of everyone's clones.

    python scripts/block_genetic_data.py <files...>
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

#: A data row from any of the consumer formats: an rsID, then a chromosome, then
#: a position, then a call. Separator varies by provider and quoting does too, so
#: this stays loose.
_SNP_ROW = re.compile(
    r"""^\s*"?(rs\d+|i\d+)"?      # rsID or 23andMe internal id
        \s*[,\t]\s*"?(\d{1,2}|[XYM]T?)"?   # chromosome
        \s*[,\t]\s*"?\d+"?        # position
        \s*[,\t]\s*"?[ACGTDI\-0]  # first allele
    """,
    re.VERBOSE | re.IGNORECASE,
)

#: Header lines that identify a provider export beyond doubt.
_PROVIDER_MARKERS = (
    "23andme",
    "ancestrydna",
    "myheritage",
    "familytreedna",
    "living dna",
    "tellmegen",
)

#: Below this, a file is a hand-written fixture rather than an export. A real
#: export has 600,000 rows; nobody types fifty by hand and calls it their genome.
_FIXTURE_LIMIT = 200

_SKIP_SUFFIXES = frozenset({".py", ".md", ".yaml", ".yml", ".toml", ".lock", ".cfg", ".ini"})


def scan(path: Path) -> tuple[int, bool]:
    """Return the number of SNP-shaped rows, and whether a provider header is present."""
    rows = 0
    provider = False

    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                lowered = line.lower()
                if line.startswith("#") and any(m in lowered for m in _PROVIDER_MARKERS):
                    provider = True
                elif _SNP_ROW.match(line):
                    rows += 1
    except OSError:
        return 0, False

    return rows, provider


def check(path: Path) -> str | None:
    if path.suffix.lower() in _SKIP_SUFFIXES or not path.is_file():
        return None

    rows, provider = scan(path)
    if rows == 0:
        return None

    # Fixtures live here and are allowed to be SNP-shaped, because that is the
    # point of them. They are still held to a size limit: a "fixture" with
    # thousands of rows is an export someone dropped in the permitted directory.
    in_fixtures = "tests/fixtures" in path.as_posix()

    if in_fixtures and rows <= _FIXTURE_LIMIT:
        return None

    if in_fixtures:
        return (
            f"{path}: {rows:,} SNP rows. Fixtures are capped at {_FIXTURE_LIMIT}.\n"
            f"    A file this large is an export, not a fixture. If a test needs "
            f"this much data, the test needs rethinking."
        )

    reason = "carries a provider header and " if provider else ""
    return (
        f"{path}: {reason}contains {rows:,} SNP-shaped rows.\n"
        f"    This looks like real genetic data. It must not be committed.\n"
        f"    See tests/fixtures/README.md."
    )


def main(argv: list[str]) -> int:
    problems = [problem for path in argv if (problem := check(Path(path))) is not None]

    if not problems:
        return 0

    print("\nBlocked: possible genetic data in this commit.\n", file=sys.stderr)
    for problem in problems:
        print(f"  {problem}\n", file=sys.stderr)
    print(
        "Genetic data cannot be revoked, and it exposes your relatives too.\n"
        "If this is genuinely a synthetic fixture, override with:\n"
        "    git commit --no-verify\n",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
