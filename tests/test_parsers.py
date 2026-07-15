"""Tests for format detection and provider parsers.

The detection tests matter more than they look. 23andMe and AncestryDNA files
are similar enough that a loose matcher will claim the wrong one, and since
AncestryDNA splits its alleles across two columns, a 23andMe parser reading an
AncestryDNA file produces a genotype of `A` where the truth is `AG`. That is a
wrong answer rather than a crash, which is the failure mode worth testing for.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from dnadistiller.models import Build
from dnadistiller.parsers import UnknownFormatError, detect_format, parse
from dnadistiller.parsers.providers import (
    AncestryDNAParser,
    FamilyTreeDNAParser,
    MyHeritageParser,
    TwentyThreeAndMeParser,
    VCFParser,
)

FIXTURES = Path(__file__).parent / "fixtures"


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------


def test_detects_23andme():
    assert detect_format(FIXTURES / "23andme_basic.txt") is TwentyThreeAndMeParser


def test_detects_ancestry():
    assert detect_format(FIXTURES / "ancestry_basic.txt") is AncestryDNAParser


def test_ancestry_is_not_claimed_by_the_23andme_parser():
    """Both are tab-separated with commented headers. Order in PARSERS matters."""
    head = (FIXTURES / "ancestry_basic.txt").read_text().splitlines()[:40]
    assert AncestryDNAParser.matches(head)
    assert not TwentyThreeAndMeParser.matches(head)


def test_unknown_format_raises_rather_than_guessing(tmp_path):
    junk = tmp_path / "junk.txt"
    junk.write_text("this is not a DNA file\njust some text\n")
    with pytest.raises(UnknownFormatError, match="Could not recognise"):
        detect_format(junk)


def test_unknown_format_error_shows_the_file_head(tmp_path):
    """The error is the bug report. It has to carry enough to act on."""
    junk = tmp_path / "junk.txt"
    junk.write_text("MYSTERY_HEADER_LINE\n")
    with pytest.raises(UnknownFormatError, match="MYSTERY_HEADER_LINE"):
        detect_format(junk)


def test_empty_file_raises(tmp_path):
    empty = tmp_path / "empty.txt"
    empty.write_text("")
    with pytest.raises(UnknownFormatError):
        detect_format(empty)


# ---------------------------------------------------------------------------
# 23andMe
# ---------------------------------------------------------------------------


def test_23andme_parses_genotypes():
    genome = parse(FIXTURES / "23andme_basic.txt")
    assert genome.source == "23andMe"
    assert genome.get("rs1801133").alleles == "AG"
    assert genome.get("rs429358").alleles == "CT"


def test_23andme_reads_build_from_header():
    assert parse(FIXTURES / "23andme_basic.txt").build is Build.GRCH37


def test_23andme_build_36_is_recognised_not_mislabelled():
    """Exports from around 2011 declare build 36, and they are still in drawers.

    Build tracks the download date, not the chip: the same v3 chip exported build
    36 in 2011 and build 37 in 2014. Mapping "anything that isn't 37" to 38 would
    report the oldest files as the newest assembly.
    """
    genome = parse(FIXTURES / "23andme_build36.txt")
    assert genome.build is Build.NCBI36


def test_build_36_positions_differ_but_rsid_lookup_still_works():
    """The reason lookups key on rsID rather than position.

    rs429358 sits at 19:50103781 in build 36 and 19:45411941 in build 37, a shift
    of about 4.7 Mb. A position-keyed tool reads the wrong locus; an rsID-keyed
    one just works, which is why build is provenance here and not a coordinate
    system we have to remap.
    """
    old = parse(FIXTURES / "23andme_build36.txt")
    new = parse(FIXTURES / "23andme_basic.txt")

    assert old.get("rs429358").position != new.get("rs429358").position
    assert old.get("rs429358") is not None
    assert new.get("rs429358") is not None


def test_build_absent_from_header_is_unknown_not_assumed():
    """Reporting a guessed build would remove the only clue an export is ancient."""
    genome = parse(FIXTURES / "23andme_five_column.txt")
    assert genome.build is Build.GRCH37  # this fixture does state it

    head = ["# rsid\tchromosome\tposition\tgenotype"]
    from dnadistiller.parsers.providers import _detect_build

    assert _detect_build(head) is Build.UNKNOWN


def test_23andme_five_column_variant_joins_alleles():
    """23andMe ships a 5-column layout too, splitting the genotype in two.

    Reading column 4 alone would turn every heterozygote into a single allele:
    `AG` reported as `A`. A wrong genotype, not a crash.
    """
    genome = parse(FIXTURES / "23andme_five_column.txt")
    assert genome.get("rs1801133").alleles == "AG"
    assert genome.get("rs429358").alleles == "CT"
    assert genome.get("rs4477212").alleles == "AA"


def test_23andme_five_column_haploid_row_has_only_one_allele():
    """A Y call in the 5-column layout has no allele2 column to join."""
    assert parse(FIXTURES / "23andme_five_column.txt").get("rs11").alleles == "A"


def test_23andme_no_call():
    assert parse(FIXTURES / "23andme_basic.txt").get("rs1799752").is_no_call


def test_23andme_internal_i_ids_are_kept():
    """23andMe assigns i-prefixed ids to variants with no rsID. They are real data."""
    assert parse(FIXTURES / "23andme_basic.txt").get("i3003137") is not None


def test_23andme_haploid_calls():
    """Y and MT report a single allele because there is only one copy."""
    genome = parse(FIXTURES / "23andme_basic.txt")
    assert genome.get("rs11").is_haploid
    assert genome.get("rs12").is_haploid


def test_23andme_indel_alleles_preserved():
    assert parse(FIXTURES / "23andme_basic.txt").get("rs113993960").alleles == "II"


# ---------------------------------------------------------------------------
# AncestryDNA
# ---------------------------------------------------------------------------


def test_ancestry_joins_split_allele_columns():
    """The bug this guards: reading allele1 only, and reporting AG as A."""
    genome = parse(FIXTURES / "ancestry_basic.txt")
    assert genome.get("rs1801133").alleles == "AG"
    assert genome.get("rs3094315").alleles == "AG"


def test_ancestry_header_row_is_not_parsed_as_data():
    genome = parse(FIXTURES / "ancestry_basic.txt")
    assert genome.get("rsid") is None


def test_ancestry_no_call_normalised():
    """AncestryDNA writes 0/0 where 23andMe writes --. Downstream sees one spelling."""
    call = parse(FIXTURES / "ancestry_basic.txt").get("rs21")
    assert call.is_no_call
    assert call.alleles == "--"


@pytest.mark.parametrize(
    ("rsid", "expected"),
    [("rs22", "X"), ("rs23", "Y"), ("rs25", "X")],
)
def test_ancestry_numeric_chromosomes_normalised(rsid, expected):
    """AncestryDNA numbers sex chromosomes PLINK-style: 23=X, 24=Y, 25=PAR, 26=MT.

    Code 26 is absent from the fixture on purpose: Ancestry V1 files carry
    chromosomes 1-25 and no mitochondrial data at all.
    """
    assert parse(FIXTURES / "ancestry_basic.txt").get(rsid).chromosome == expected


def test_ancestry_crlf_homozygote_is_not_read_as_heterozygote():
    """Real AncestryDNA exports use CRLF, and this is the trap in them.

    With a naive line split, allele2 keeps its trailing carriage return: `A` and
    `A\\r` are unequal, so every homozygote silently reads as a heterozygote.
    Nothing errors. The zygosity is just wrong, for every marker, in the
    direction that invents findings.
    """
    genome = parse(FIXTURES / "ancestry_crlf.txt")

    call = genome.get("rs4477212")
    assert call.alleles == "AA"
    assert "\r" not in call.alleles
    assert call.canonical() == "AA"


def test_ancestry_crlf_does_not_corrupt_the_last_column():
    genome = parse(FIXTURES / "ancestry_crlf.txt")
    for call in genome:
        assert "\r" not in call.alleles, f"{call.rsid} kept a carriage return"


def test_ancestry_reversed_allele_order_canonicalises():
    """Real Ancestry files contain the same het call written both ways.

    One 2013 export had 44,728 `AG` and 39,678 `GA`. Order carries no meaning:
    an array does not resolve which chromosome an allele sits on.
    """
    genome = parse(FIXTURES / "ancestry_crlf.txt")
    assert genome.get("rs3131972").canonical() == genome.get("rs3094315").canonical() == "AG"


def test_apoe_is_not_computable_from_ancestry():
    """AncestryDNA V1 carries rs7412 but not rs429358.

    This is the real "missing APOE" case, and it belongs to Ancestry rather than
    23andMe. Half an APOE result is worse than none, so the derivation must
    decline rather than reason from rs7412 alone.
    """
    from dnadistiller.derive import derive_apoe
    from dnadistiller.models import FindingStatus

    from .test_derive import apoe_marker

    genome = parse(FIXTURES / "ancestry_basic.txt")
    assert genome.get("rs7412") is not None
    assert genome.get("rs429358") is None

    assert derive_apoe(genome, apoe_marker()).status is FindingStatus.NOT_ON_CHIP


# ---------------------------------------------------------------------------
# CSV providers
# ---------------------------------------------------------------------------


def test_myheritage(tmp_path):
    path = tmp_path / "mh.csv"
    path.write_text(
        "# MyHeritage DNA raw data.\n"
        "# For more information visit: https://www.myheritage.com/dna\n"
        "RSID,CHROMOSOME,POSITION,RESULT\n"
        '"rs1801133","1","11856378","AG"\n'
        '"rs7412","19","45412079","CC"\n'
    )
    assert detect_format(path) is MyHeritageParser
    genome = parse(path)
    assert genome.get("rs1801133").alleles == "AG"
    assert genome.source == "MyHeritage"


def test_ftdna_without_comment_header(tmp_path):
    """FTDNA ships a bare quoted CSV whose first line is the column header."""
    path = tmp_path / "ftdna.csv"
    path.write_text('RSID,CHROMOSOME,POSITION,RESULT\n"rs1801133","1","11856378","AG"\n')
    assert detect_format(path) is FamilyTreeDNAParser
    assert parse(path).get("rs1801133").alleles == "AG"


# ---------------------------------------------------------------------------
# VCF
# ---------------------------------------------------------------------------


def vcf(body: str) -> str:
    return (
        "##fileformat=VCFv4.2\n"
        "##reference=GRCh37\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n" + body
    )


def test_vcf_heterozygous(tmp_path):
    path = tmp_path / "s.vcf"
    path.write_text(vcf("1\t11856378\trs1801133\tG\tA\t.\tPASS\t.\tGT\t0/1\n"))
    assert detect_format(path) is VCFParser
    assert parse(path).get("rs1801133").alleles == "AG"


def test_vcf_homozygous_alt(tmp_path):
    path = tmp_path / "s.vcf"
    path.write_text(vcf("1\t11856378\trs1801133\tG\tA\t.\tPASS\t.\tGT\t1/1\n"))
    assert parse(path).get("rs1801133").alleles == "AA"


def test_vcf_phased_separator(tmp_path):
    path = tmp_path / "s.vcf"
    path.write_text(vcf("1\t11856378\trs1801133\tG\tA\t.\tPASS\t.\tGT\t0|1\n"))
    assert parse(path).get("rs1801133").alleles == "AG"


def test_vcf_skips_indels_rather_than_misreading_them(tmp_path):
    """This reader handles SNVs. An indel it cannot read must be absent, not wrong."""
    path = tmp_path / "s.vcf"
    path.write_text(vcf("1\t100\trs999\tGA\tG\t.\tPASS\t.\tGT\t0/1\n"))
    assert parse(path).get("rs999") is None


def test_vcf_skips_missing_calls(tmp_path):
    path = tmp_path / "s.vcf"
    path.write_text(vcf("1\t100\trs999\tG\tA\t.\tPASS\t.\tGT\t./.\n"))
    assert parse(path).get("rs999") is None


def test_vcf_skips_rows_without_an_rsid(tmp_path):
    path = tmp_path / "s.vcf"
    path.write_text(vcf("1\t100\t.\tG\tA\t.\tPASS\t.\tGT\t0/1\n"))
    assert len(parse(path)) == 0


# ---------------------------------------------------------------------------
# Regressions drawn from failure modes observed in another implementation.
#
# Each of these is a real bug found in a shipped consumer-DNA tool. They are
# worth pinning here precisely because none of them raise: each one produces a
# confident, wrong, silent result on a file the tool claims to support.
# ---------------------------------------------------------------------------


def test_space_separated_file_still_parses():
    """Headers say TAB-separated; files in the wild sometimes are not.

    Exports exist whose header insists "Fields are TAB-separated" while the data
    below is space-separated, mangled by some tool in between. A tab-only split
    does not degrade on these, it fails totally: every row becomes one field.
    """
    from dnadistiller.parsers.providers import _split

    assert _split("rs1\t1\t82154\tAA") == ["rs1", "1", "82154", "AA"]
    assert _split("rs1 1 82154 AA") == ["rs1", "1", "82154", "AA"]


def test_tab_split_wins_over_whitespace_fallback():
    """The fallback must not corrupt a genuine tab row containing a space."""
    from dnadistiller.parsers.providers import _split

    assert _split("rs1\t1\t82154\tA A") == ["rs1", "1", "82154", "A A"]


def test_i_ids_are_not_discarded_as_malformed():
    """23andMe assigns i-prefixed ids to variants with no rsID.

    A real v3 file has around 10,000 of them, and they include CFTR and BRCA
    markers. Filtering rows on `rsid.startswith("rs")` throws all of them away
    and counts them as parse errors.
    """
    genome = parse(FIXTURES / "23andme_basic.txt")
    assert genome.get("i3003137") is not None
    assert genome.malformed_rows == 0


def test_ancestry_five_columns_are_not_read_as_four():
    """Reading fields[3] as the genotype on a 5-column file halves every het.

    `A` instead of `AG`. No error, no warning, wrong zygosity on every marker.
    """
    genome = parse(FIXTURES / "ancestry_basic.txt")
    assert genome.get("rs3094315").alleles == "AG"
    assert len(genome.get("rs3094315").alleles) == 2


def test_bom_does_not_break_the_comment_header():
    """A UTF-8 BOM glues itself to the first character of line 1.

    With plain utf-8 the header comment no longer starts with '#', so it parses
    as data, and the build declaration it carries is lost.
    """
    import tempfile

    content = (
        "# This data file generated by 23andMe at: x\n"
        "# reference human assembly build 37\n"
        "# rsid\tchromosome\tposition\tgenotype\n"
        "rs1801133\t1\t11856378\tAG\n"
    )
    with tempfile.NamedTemporaryFile("wb", suffix=".txt", delete=False) as handle:
        handle.write(b"\xef\xbb\xbf" + content.encode())
        path = Path(handle.name)

    genome = parse(path, expect_full_export=False)
    assert genome.build is Build.GRCH37, "BOM swallowed the build declaration"
    assert genome.get("rs1801133").alleles == "AG"
    assert genome.malformed_rows == 0
    path.unlink()
