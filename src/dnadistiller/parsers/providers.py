"""One parser per provider export format.

Adding a provider means adding a class here and listing it in `PARSERS`. Keep
`matches` narrow: these formats resemble each other closely enough that a loose
matcher will happily claim a competitor's file and misread every column in it.
"""

from __future__ import annotations

import csv
import re
from abc import ABC, abstractmethod
from pathlib import Path

from ..models import NO_CALL_TOKENS, Build, Genome, Genotype, ParseIssue, Severity

#: 23andMe and AncestryDNA both state their assembly in the header comments,
#: e.g. "reference human assembly build 37". Always read it, never assume it:
#: 23andMe exports from around 2011 say build 36, and the same v3 chip produced
#: build 36 in 2011 and build 37 in 2014. Build tracks the download date, not the
#: chip version, so there is nothing else to infer it from.
_BUILD_RE = re.compile(r"build\s+(3[678])\b", re.IGNORECASE)

_BUILDS: dict[str, Build] = {
    "36": Build.NCBI36,
    "37": Build.GRCH37,
    "38": Build.GRCH38,
}

#: AncestryDNA numbers the sex chromosomes PLINK-style instead of naming them.
#: 25 is the pseudoautosomal region of X, which is diploid in everyone and so is
#: reported separately from X proper.
_NUMERIC_CHROMOSOMES: dict[str, str] = {
    "23": "X",
    "24": "Y",
    "25": "X",
    "26": "MT",
}


def _split(line: str) -> list[str]:
    """Split a data row on tabs, falling back to arbitrary whitespace.

    Every tab-separated provider says TAB-separated in its own header and means
    it. Files still arrive with the tabs converted to spaces, because they have
    passed through a spreadsheet, a mail client, or a copy-paste. Exports exist
    in the wild whose header insists "Fields are TAB-separated" while the data
    below is space-separated.

    The fallback only engages when the line contains no tab at all, so it cannot
    corrupt a genuine tab-separated row that happens to hold a space. Without it
    such a file does not degrade, it fails completely: every row splits to one
    field and the whole export reads as malformed.
    """
    stripped = line.rstrip("\r\n")
    if "\t" in stripped:
        return stripped.split("\t")
    return stripped.split()


def _normalise_chromosome(raw: str) -> str:
    """Bring chromosome labels to one convention: 1-22, X, Y, MT.

    Providers variously write MT, M, chrM, or 26 for the mitochondrion. Since
    findings are keyed on rsID rather than locus this is mostly cosmetic — but
    it is what the `parse` sanity check displays, and a distribution split
    across `MT` and `26` looks like a bug to the person reading it.
    """
    cleaned = raw.strip().upper().removeprefix("CHR")
    if cleaned in _NUMERIC_CHROMOSOMES:
        return _NUMERIC_CHROMOSOMES[cleaned]
    if cleaned in {"M", "MT"}:
        return "MT"
    return cleaned


def _detect_build(head: list[str]) -> Build:
    """Read the assembly out of the header comments.

    Returns UNKNOWN rather than guessing 37 when the header does not say. A
    wrong build here is only provenance being misreported, since lookups key on
    rsID, but misreporting it would undermine the one signal a user has that
    their export is ancient.
    """
    for line in head:
        match = _BUILD_RE.search(line)
        if match:
            return _BUILDS[match.group(1)]
    return Build.UNKNOWN


#: Real exports carry 500,000 to 700,000 variants. Well under that means the file
#: is truncated, filtered, or a fragment someone pasted together, and every
#: "not tested" it produces is an artefact rather than a fact about the chip.
#: Silent success on such a file is the worst outcome available: it looks exactly
#: like a clean run that happened to find nothing.
_EXPECTED_MINIMUM_VARIANTS = 100_000

#: Above this, a parse is not merely imperfect, it is wrong. Malformed rows in
#: the low hundreds are normal; a large fraction means the delimiter or the
#: column layout was misread, and the calls that *did* parse are not trustworthy
#: either.
_MALFORMED_ROW_ALARM = 0.02

#: Ambiguous calls cluster when an encoding or delimiter problem is mangling the
#: last column. Biology does not produce them at this rate.
_AMBIGUOUS_CALL_ALARM = 0.02


def diagnose(genome: Genome, *, expect_full_export: bool = True) -> list[ParseIssue]:
    """Judge whether a parse should be trusted, and say so in structured form.

    Separated from parsing because the parsers are per-provider and these
    judgements are not: a truncated file is truncated whoever made it. Returns
    issues rather than raising, since none of these are fatal and all of them are
    the user's call.

    `expect_full_export` is off for fixtures and fragments, where a small file is
    the point rather than a symptom.
    """
    issues: list[ParseIssue] = []
    total = len(genome)

    if total == 0:
        issues.append(
            ParseIssue(
                Severity.ERROR,
                "no-variants",
                "The file matched a known format but produced no variants. It is "
                "probably truncated or the columns are not where we expect them.",
            )
        )
        return issues

    if expect_full_export and total < _EXPECTED_MINIMUM_VARIANTS:
        issues.append(
            ParseIssue(
                Severity.WARNING,
                "small-file",
                f"Only {total:,} variants. A real export carries 500,000 or more, so "
                "this file looks truncated or filtered. Treat every 'not tested' "
                "below as a fact about the file rather than about your chip.",
            )
        )

    rows = total + genome.malformed_rows
    if genome.malformed_rows and rows and genome.malformed_rows / rows > _MALFORMED_ROW_ALARM:
        issues.append(
            ParseIssue(
                Severity.ERROR,
                "malformed-rows",
                f"Skipped {genome.malformed_rows:,} unreadable rows out of {rows:,}. "
                "At this rate the format was probably misdetected, which means the "
                "rows that did parse may be misread too. Please report this file's "
                "first few lines.",
            )
        )
    elif genome.malformed_rows:
        issues.append(
            ParseIssue(
                Severity.INFO,
                "malformed-rows",
                f"Skipped {genome.malformed_rows:,} unreadable rows. A few is normal.",
            )
        )

    if genome.duplicate_rsids:
        issues.append(
            ParseIssue(
                Severity.INFO,
                "duplicate-rsids",
                f"{genome.duplicate_rsids:,} rsIDs appeared more than once. The last "
                "occurrence was kept. This is common in real exports.",
            )
        )

    if genome.ambiguous_calls and genome.ambiguous_calls / total > _AMBIGUOUS_CALL_ALARM:
        issues.append(
            ParseIssue(
                Severity.WARNING,
                "ambiguous-calls",
                f"{genome.ambiguous_calls:,} calls are in a format we do not recognise "
                "as either a genotype or a no-call. A rate this high usually means an "
                "encoding or delimiter problem rather than unusual biology.",
            )
        )

    if genome.source == GenericParser.provider:
        issues.append(
            ParseIssue(
                Severity.WARNING,
                "unknown-vendor",
                "This file's vendor could not be identified. It was read with a generic "
                "layout that matched its structure, so the genotypes are as the file "
                "states them, but vendor-specific corrections were not applied. Please "
                "report the file's header lines so the vendor can be supported properly.",
            )
        )

    if genome.build is Build.UNKNOWN:
        issues.append(
            ParseIssue(
                Severity.INFO,
                "build-unknown",
                "The file does not state its reference build. Harmless here, because "
                "markers are found by rsID rather than by position.",
            )
        )

    if genome.build is Build.NCBI36:
        issues.append(
            ParseIssue(
                Severity.INFO,
                "old-build",
                "This export uses NCBI36, which providers stopped issuing around 2013. "
                "It still works, but a fresh download would cover more of the catalog.",
            )
        )

    no_calls = sum(1 for g in genome if g.is_no_call)
    rate = no_calls / total
    if rate > 0.05:
        issues.append(
            ParseIssue(
                Severity.WARNING,
                "high-no-call-rate",
                f"{rate:.1%} of calls failed. Above about 5% the sample quality is "
                "questionable. Consider re-downloading, or asking the provider to "
                "reprocess the sample.",
            )
        )

    return issues


#: What every consumer array structurally cannot see. Attached to every parse so
#: the limitation travels with the results rather than living in a FAQ nobody
#: reads. The point these make is the one users most reliably get wrong: a normal
#: result at a tested position says nothing about the untested rest of the gene.
GENERIC_BLIND_SPOTS: tuple[str, ...] = (
    "A consumer array assays a tiny, preselected fraction of the genome. A normal "
    "result here never rules out other variants in the same gene.",
    "This file type cannot resolve rare variants, structural variants, copy-number "
    "changes, methylation, HLA types, or most CYP2D6 star alleles. Those need "
    "sequencing.",
)


class Parser(ABC):
    """Base class for provider parsers."""

    #: Human-readable name, shown in error messages and in profile provenance.
    provider: str = "unknown"

    #: Limitations specific to this provider, added to GENERIC_BLIND_SPOTS.
    #: Kept per-parser rather than in a central table so that adding a provider
    #: forces an answer to "what can this one not see?" at the same moment as
    #: "how do I read it?".
    blind_spots: tuple[str, ...] = ()

    @classmethod
    @abstractmethod
    def matches(cls, head: list[str]) -> bool:
        """Whether this parser recognises a file from its first lines."""

    @abstractmethod
    def parse(self, path: Path) -> Genome:
        """Read the file into a Genome."""

    def _new_genome(self) -> Genome:
        return Genome(
            source=self.provider, blind_spots=list(GENERIC_BLIND_SPOTS + self.blind_spots)
        )


class TwentyThreeAndMeParser(Parser):
    """23andMe raw export (v3, v4, v5).

    Tab-separated, with `#` comment lines carrying the metadata. The final
    comment line is the column header, which is why skipping comments and
    supplying our own column names works rather than being a shortcut::

        # This data file generated by 23andMe at: ...
        # reference human assembly build 37
        # rsid	chromosome	position	genotype
        rs4477212	1	82154	AA

    Calls on a male X or Y outside the pseudoautosomal regions, and on MT,
    arrive as a single letter because there is only one copy to report.
    """

    provider = "23andMe"
    blind_spots = (
        "23andMe chip content changed between v3, v4, and v5. A missing marker here "
        "usually reflects which chip you were tested on rather than a real negative.",
        "Sample-swap and no-call rates are low but not zero. A single surprising result "
        "is worth confirming before you believe it.",
    )

    @classmethod
    def matches(cls, head: list[str]) -> bool:
        text = "\n".join(head)
        if "23andme" not in text.lower():
            return False

        # Require the column header, and require it to be commented. Both facts
        # do work: AncestryDNA files mention 23andMe in prose about imports, and
        # the 5-column 23andMe layout has the same allele1/allele2 columns as
        # AncestryDNA. What separates them is that 23andMe comments its header
        # and AncestryDNA does not.
        return any(
            line.startswith("#")
            and "rsid" in line.lower()
            and ("genotype" in line.lower() or "allele1" in line.lower())
            for line in head
        )

    def parse(self, path: Path) -> Genome:
        genome = self._new_genome()
        head: list[str] = []

        with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    if len(head) < 40:
                        head.append(line)
                    continue

                fields = _split(line)
                if len(fields) < 4:
                    genome.malformed_rows += 1
                    continue

                rsid = fields[0].strip()
                if rsid.lower() == "rsid":  # uncommented header, seen in some exports
                    continue

                try:
                    pos = int(fields[2])
                except ValueError:
                    genome.malformed_rows += 1
                    continue

                # 23andMe ships two layouts. The usual one puts the genotype in a
                # single column; a 5-column variant splits it across allele1 and
                # allele2 the way AncestryDNA does. Taking fields[3] blindly would
                # read `AG` as `A` on those files: a wrong genotype, silently, with
                # every heterozygote reduced to whichever allele came first.
                alleles = "".join(f.strip() for f in fields[3:5]) if len(fields) >= 5 else fields[3]

                genome.add(
                    Genotype(
                        rsid=rsid,
                        chromosome=_normalise_chromosome(fields[1]),
                        position=pos,
                        alleles=alleles.strip().upper(),
                    )
                )

        genome.build = _detect_build(head)
        return genome


class AncestryDNAParser(Parser):
    """AncestryDNA raw export.

    Differs from 23andMe in two ways that matter: an uncommented column-header
    row, and alleles split across two columns which we rejoin. A no-call is `0`
    in each allele column rather than `--`::

        #AncestryDNA raw data download
        rsid	chromosome	position	allele1	allele2
        rs4477212	1	82154	A	A
    """

    provider = "AncestryDNA"
    blind_spots = (
        "AncestryDNA does not genotype rs429358, so APOE cannot be determined from this "
        "file at all. That is a gap in the chip, not a reassuring result.",
        "AncestryDNA V1 files carry no mitochondrial data.",
        "Pharmacogenomic coverage is patchier than 23andMe's. A missing drug marker is "
        "chip design rather than biology.",
    )

    @classmethod
    def matches(cls, head: list[str]) -> bool:
        text = "\n".join(head).lower()
        if "ancestrydna" in text:
            return True
        # Fall back to the structure: split allele columns are theirs alone.
        return any(
            "allele1" in line.lower() and "allele2" in line.lower()
            for line in head
            if not line.startswith("#")
        )

    def parse(self, path: Path) -> Genome:
        genome = self._new_genome()
        head: list[str] = []
        seen_header = False

        with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    if len(head) < 40:
                        head.append(line)
                    continue

                fields = _split(line)
                if not seen_header:
                    # The first uncommented line is the column header.
                    seen_header = True
                    if fields and fields[0].strip().lower() == "rsid":
                        continue

                if len(fields) < 5:
                    genome.malformed_rows += 1
                    continue

                rsid, chromosome, position, allele1, allele2 = fields[:5]
                try:
                    pos = int(position)
                except ValueError:
                    genome.malformed_rows += 1
                    continue

                alleles = f"{allele1.strip()}{allele2.strip()}".upper()
                # "00" is AncestryDNA's no-call; rewrite to the common "--" so
                # downstream code has one spelling to recognise.
                if alleles == "00":
                    alleles = "--"

                genome.add(
                    Genotype(
                        rsid=rsid.strip(),
                        chromosome=_normalise_chromosome(chromosome),
                        position=pos,
                        alleles=alleles,
                    )
                )

        genome.build = _detect_build(head)
        return genome


class _QuotedCSVParser(Parser):
    """Shared implementation for the CSV-with-quoted-fields providers.

    MyHeritage and FamilyTreeDNA emit effectively the same file — header
    `RSID,CHROMOSOME,POSITION,RESULT`, every field quoted — and differ only in
    whether comments precede it. They stay separate classes so that provenance
    in the profile names the actual provider.
    """

    provider = "csv"

    def parse(self, path: Path) -> Genome:
        genome = self._new_genome()
        head: list[str] = []

        with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
            rows = (line for line in handle if not self._is_comment(line, head))
            for fields in csv.reader(rows):
                if len(fields) < 4:
                    continue
                rsid, chromosome, position, result = (f.strip() for f in fields[:4])
                if rsid.lower() == "rsid":
                    continue
                try:
                    pos = int(position)
                except ValueError:
                    continue

                genome.add(
                    Genotype(
                        rsid=rsid,
                        chromosome=_normalise_chromosome(chromosome),
                        position=pos,
                        alleles=result.upper(),
                    )
                )

        genome.build = _detect_build(head)
        return genome

    @staticmethod
    def _is_comment(line: str, head: list[str]) -> bool:
        if line.startswith("#"):
            if len(head) < 40:
                head.append(line)
            return True
        return False


class MyHeritageParser(_QuotedCSVParser):
    provider = "MyHeritage"
    blind_spots = (
        "MyHeritage coverage is strongest for common ancestry-informative SNPs. Niche "
        "pharmacogenomic markers are patchier.",
    )

    @classmethod
    def matches(cls, head: list[str]) -> bool:
        return "myheritage" in "\n".join(head).lower()


class FamilyTreeDNAParser(_QuotedCSVParser):
    provider = "FamilyTreeDNA"
    blind_spots = (
        "FamilyTreeDNA's autosomal chip is built for genealogy. A missing rsID here is "
        "vendor omission rather than a negative genotype.",
    )

    #: MyHeritage and FTDNA both use this header. MyHeritage V1 writes RESULT and
    #: V2 writes RESULTS, so match the stem rather than the exact string: an
    #: exact match silently stopped recognising V2 files the day it shipped.
    _HEADER_STEM = "RSID,CHROMOSOME,POSITION,RESULT"

    @classmethod
    def matches(cls, head: list[str]) -> bool:
        text = "\n".join(head).lower()
        if "familytreedna" in text or "ftdna" in text:
            return True
        # FTDNA ships no comment header at all: an unadorned quoted CSV whose
        # first line is the column header. Match that shape, but only when no
        # other provider has claimed the file — hence its late position in PARSERS.
        first = next((line for line in head if line.strip()), "")
        return first.replace('"', "").strip().upper().startswith(cls._HEADER_STEM)


class LivingDNAParser(Parser):
    """Living DNA raw export: tab-separated, 23andMe-shaped, own header wording."""

    provider = "LivingDNA"

    @classmethod
    def matches(cls, head: list[str]) -> bool:
        return "living dna" in "\n".join(head).lower() or "livingdna" in "\n".join(head).lower()

    def parse(self, path: Path) -> Genome:
        genome = TwentyThreeAndMeParser().parse(path)
        genome.source = self.provider
        return genome


class GenericParser(Parser):
    """Last resort for an export whose vendor we do not recognise.

    Consumer genotyping is a long tail. Beyond the handful of vendors named
    above there are tellmeGen, Genes for Good, Living DNA rebrands, regional
    labs, and whatever launched last month, and they nearly all emit the same
    shape: an identifier, a chromosome, a position, and a call. Refusing every
    one of them because the header does not say a name we know serves nobody.

    This parser is the exception to "refuse rather than guess", and it is
    allowed to exist only because it does not guess. The rule it breaks was
    never really about vendors: it was about not misreading columns. So instead
    of trusting a vendor name, `matches` proves the layout from the data itself,
    over many rows:

      * column 0 looks like an rsID or a 23andMe i-id
      * column 2 parses as an integer
      * the final column is allele-shaped, or a recognised no-call

    A file where all three hold across most of its rows has its columns where we
    think they are, whoever made it. A file where they do not is rejected, and
    the user gets the unknown-format error as before.

    It is still last in PARSERS: any named vendor claims its own file first,
    because a specific parser knows things about its format that shape alone
    cannot show, such as AncestryDNA's PLINK chromosome coding.
    """

    provider = "unknown vendor"
    blind_spots = (
        "This file's vendor could not be identified, so it was read using a generic "
        "layout that matched its structure. The genotypes are as the file states "
        "them, but vendor-specific quirks cannot be corrected for. Treat the results "
        "as provisional and please report the file's header so the vendor can be "
        "supported properly.",
    )

    #: Fraction of sampled data rows that must fit the expected shape. Set high:
    #: a partial match means the layout is not what we think, and a wrong column
    #: guess produces confident wrong genotypes rather than an error.
    _SHAPE_THRESHOLD = 0.9

    #: Minimum data rows needed before the shape argument means anything. Three
    #: rows can fit almost any layout by chance.
    _MIN_ROWS = 5

    @classmethod
    def matches(cls, head: list[str]) -> bool:
        rows = [line for line in head if line.strip() and not line.startswith("#")]
        if len(rows) < cls._MIN_ROWS:
            return False

        fits = sum(1 for row in rows if cls._row_has_expected_shape(_split(row)))
        return fits / len(rows) >= cls._SHAPE_THRESHOLD

    @staticmethod
    def _row_has_expected_shape(fields: list[str]) -> bool:
        if len(fields) < 4:
            return False

        rsid = fields[0].strip().lower()
        if not (rsid.startswith("rs") and rsid[2:].isdigit()) and not (
            rsid.startswith("i") and rsid[1:].isdigit()
        ):
            return False

        if not fields[2].strip().isdigit():
            return False

        call = "".join(f.strip() for f in fields[3:5]).upper()
        if call in NO_CALL_TOKENS or set(call) <= {"-", "0"}:
            return True
        return 1 <= len(call) <= 2 and all(allele in "ACGTID" for allele in call)

    def parse(self, path: Path) -> Genome:
        genome = self._new_genome()
        head: list[str] = []

        with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    if len(head) < 40:
                        head.append(line)
                    continue

                fields = _split(line)
                if len(fields) < 4:
                    genome.malformed_rows += 1
                    continue

                rsid = fields[0].strip()
                if rsid.lower() == "rsid":  # a column header, whoever wrote it
                    continue

                try:
                    pos = int(fields[2])
                except ValueError:
                    genome.malformed_rows += 1
                    continue

                # Same two layouts everyone else uses: one genotype column, or a
                # split allele1/allele2 pair.
                alleles = "".join(f.strip() for f in fields[3:5]) if len(fields) >= 5 else fields[3]
                if alleles.upper() in {"00", "0"}:
                    alleles = "--"

                genome.add(
                    Genotype(
                        rsid=rsid,
                        chromosome=_normalise_chromosome(fields[1]),
                        position=pos,
                        alleles=alleles.strip().upper(),
                    )
                )

        genome.build = _detect_build(head)
        return genome


class VCFParser(Parser):
    """Minimal VCF reader for providers that offer a VCF export.

    Deliberately narrow: single-sample, diploid, SNVs only. VCF is a large
    format and this is not a general implementation — it reads enough to look up
    catalog markers and skips anything it cannot read unambiguously, on the
    principle that a skipped marker reports "not tested" while a
    misinterpreted one reports a wrong genotype.
    """

    provider = "VCF"

    @classmethod
    def matches(cls, head: list[str]) -> bool:
        return any(line.startswith("##fileformat=VCF") for line in head)

    def parse(self, path: Path) -> Genome:
        genome = self._new_genome()
        head: list[str] = []

        with path.open("r", encoding="utf-8-sig", errors="replace") as handle:
            for line in handle:
                if line.startswith("#"):
                    if len(head) < 40:
                        head.append(line)
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue

                chromosome, position, rsid, ref, alt = fields[:5]
                sample = fields[9]

                if rsid == "." or not rsid.lower().startswith("rs"):
                    continue

                alleles = self._genotype_from_sample(sample, ref, alt)
                if alleles is None:
                    continue

                try:
                    pos = int(position)
                except ValueError:
                    continue

                genome.add(
                    Genotype(
                        rsid=rsid.strip(),
                        chromosome=_normalise_chromosome(chromosome),
                        position=pos,
                        alleles=alleles,
                    )
                )

        genome.build = _detect_build(head)
        return genome

    @staticmethod
    def _genotype_from_sample(sample: str, ref: str, alt: str) -> str | None:
        """Resolve a GT field like `0/1` into allele letters.

        Returns None for anything not a clean diploid or haploid SNV call —
        indels, multi-allelic sites, and missing calls all land here.
        """
        gt = sample.split(":")[0].replace("|", "/")
        if not gt or gt in {"./.", "."}:
            return None

        options = [ref, *alt.split(",")]
        if any(len(option) != 1 for option in options):
            return None  # indel or structural variant

        letters = []
        for index in gt.split("/"):
            if not index.isdigit() or int(index) >= len(options):
                return None
            letters.append(options[int(index)])

        return "".join(sorted(letters)).upper()
