"""Tests for parse diagnostics and vendor blind spots.

These exist because the dangerous parser failure is not the loud one. A file
that fails to open gets reported by the user. A file that parses to 4,000
variants instead of 600,000, or whose delimiter was misread so every row is
skipped, produces a profile full of "not tested" that looks exactly like a
clean run finding nothing notable.

The whole point of this module is to make that state visible.
"""

from __future__ import annotations

import pytest

from dnadistiller.models import Build, Genome, Genotype, Severity
from dnadistiller.parsers import parse
from dnadistiller.parsers.providers import diagnose


def codes(genome: Genome, *, expect_full_export: bool = True) -> set[str]:
    return {i.code for i in diagnose(genome, expect_full_export=expect_full_export)}


def genome_with(count: int, *, alleles: str = "AG") -> Genome:
    genome = Genome(source="test", build=Build.GRCH37)
    for i in range(count):
        genome.add(Genotype(rsid=f"rs{i}", chromosome="1", position=i, alleles=alleles))
    return genome


# ---------------------------------------------------------------------------
# Truncated files
# ---------------------------------------------------------------------------


def test_small_file_is_flagged():
    """The failure that looks like success.

    A truncated export yields "not tested" for most of the catalog, which reads
    as an all-clear rather than as a broken input.
    """
    assert "small-file" in codes(genome_with(500))


def test_full_size_export_is_not_flagged():
    assert "small-file" not in codes(genome_with(150_000))


def test_small_file_check_can_be_suppressed_for_fragments():
    assert "small-file" not in codes(genome_with(500), expect_full_export=False)


def test_empty_genome_is_an_error_not_a_warning():
    issues = diagnose(Genome(source="test"))
    assert [i.code for i in issues] == ["no-variants"]
    assert issues[0].severity is Severity.ERROR


# ---------------------------------------------------------------------------
# Malformed rows: the delimiter-misread signature
# ---------------------------------------------------------------------------


def test_a_few_malformed_rows_is_only_informational():
    genome = genome_with(150_000)
    genome.malformed_rows = 30
    issues = {i.code: i for i in diagnose(genome)}
    assert issues["malformed-rows"].severity is Severity.INFO


def test_mass_malformed_rows_is_an_error():
    """A high skip rate means the format was misdetected.

    That matters more than the skipped rows themselves: if the columns were
    misread, the rows that *did* parse are misread too, and those become
    findings.
    """
    genome = genome_with(150_000)
    genome.malformed_rows = 50_000
    issues = {i.code: i for i in diagnose(genome)}
    assert issues["malformed-rows"].severity is Severity.ERROR


# ---------------------------------------------------------------------------
# Duplicate and ambiguous calls
# ---------------------------------------------------------------------------


def test_duplicate_rsids_are_counted_not_hidden():
    """Last-one-wins is arbitrary, so it gets reported rather than assumed."""
    genome = Genome(source="test")
    genome.add(Genotype(rsid="rs1", chromosome="1", position=1, alleles="AA"))
    genome.add(Genotype(rsid="rs1", chromosome="1", position=1, alleles="AG"))

    assert genome.duplicate_rsids == 1
    assert genome.get("rs1").alleles == "AG"  # last wins
    assert "duplicate-rsids" in codes(genome, expect_full_export=False)


def test_duplicate_detection_is_case_insensitive():
    genome = Genome(source="test")
    genome.add(Genotype(rsid="rs1", chromosome="1", position=1, alleles="AA"))
    genome.add(Genotype(rsid="RS1", chromosome="1", position=1, alleles="AG"))
    assert genome.duplicate_rsids == 1


@pytest.mark.parametrize("alleles", ["AG", "A", "--", "II", "DD"])
def test_legitimate_calls_are_not_ambiguous(alleles):
    genotype = Genotype(rsid="rs1", chromosome="1", position=1, alleles=alleles)
    assert not genotype.is_ambiguous


@pytest.mark.parametrize("alleles", ["A\r", "AGG", "XY", "A G"])
def test_mangled_calls_are_ambiguous(alleles):
    """The signature of an encoding or delimiter problem, not of biology."""
    genotype = Genotype(rsid="rs1", chromosome="1", position=1, alleles=alleles)
    assert genotype.is_ambiguous


def test_mass_ambiguous_calls_are_flagged():
    genome = genome_with(150_000, alleles="XY")
    assert "ambiguous-calls" in codes(genome)


def test_a_stray_ambiguous_call_is_not_flagged():
    genome = genome_with(150_000)
    genome.add(Genotype(rsid="rs_odd", chromosome="1", position=1, alleles="XY"))
    assert "ambiguous-calls" not in codes(genome)


# ---------------------------------------------------------------------------
# Build reporting
# ---------------------------------------------------------------------------


def test_unknown_build_is_informational_not_alarming():
    """Harmless, because lookups key on rsID. Worth saying, not worth shouting."""
    genome = genome_with(150_000)
    genome.build = Build.UNKNOWN
    issues = {i.code: i for i in diagnose(genome)}
    assert issues["build-unknown"].severity is Severity.INFO


def test_old_build_is_reported():
    genome = genome_with(150_000)
    genome.build = Build.NCBI36
    assert "old-build" in codes(genome)


def test_high_no_call_rate_is_flagged():
    genome = genome_with(150_000, alleles="--")
    assert "high-no-call-rate" in codes(genome)


# ---------------------------------------------------------------------------
# Blind spots
# ---------------------------------------------------------------------------


def test_every_parse_carries_generic_blind_spots(tmp_path):
    genome = parse_fixture("23andme_basic.txt")
    joined = " ".join(genome.blind_spots).lower()
    assert "preselected fraction" in joined
    assert "never rules out other variants in the same gene" in joined


def test_ancestry_blind_spots_name_the_apoe_gap():
    """The most consequential vendor-specific gap in the catalog.

    An AncestryDNA user asking about APOE gets "not tested", which without this
    note is indistinguishable from a reassuring result.
    """
    genome = parse_fixture("ancestry_basic.txt")
    joined = " ".join(genome.blind_spots)
    assert "rs429358" in joined
    assert "not a reassuring result" in joined


def test_blind_spots_are_vendor_specific():
    ancestry = " ".join(parse_fixture("ancestry_basic.txt").blind_spots)
    twentythree = " ".join(parse_fixture("23andme_basic.txt").blind_spots)

    assert "rs429358" in ancestry
    assert "rs429358" not in twentythree
    assert "v3, v4, and v5" in twentythree


def parse_fixture(name: str) -> Genome:
    from pathlib import Path

    return parse(Path(__file__).parent / "fixtures" / name)
