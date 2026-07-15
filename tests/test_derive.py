"""Tests for APOE haplotype derivation.

APOE is the marker with the largest effect in the catalog and the one people
care most about, so a wrong call here matters more than anywhere else. The
diplotype table is enumerated exhaustively rather than sampled.
"""

from __future__ import annotations

import pytest

from dnadistiller.derive import derive_apoe
from dnadistiller.interpret import interpret_marker
from dnadistiller.models import (
    Category,
    Evidence,
    FindingStatus,
    Genome,
    Genotype,
    Interpretation,
    Marker,
    Zygosity,
)

ALL_DIPLOTYPES = [
    "e1/e1",
    "e1/e2",
    "e1/e4",
    "e2/e2",
    "e2/e3",
    "e2/e4",
    "e3/e3",
    "e3/e4",
    "e4/e4",
]


def apoe_marker() -> Marker:
    return Marker(
        id="apoe",
        gene="APOE",
        name="APOE isoform",
        category=Category.NEURO,
        evidence=Evidence.STRONG,
        rsids=("rs429358", "rs7412"),
        citations=("PMID:1",),
        derivation="apoe",
        sensitive=True,
        interpretations={d: Interpretation(summary=f"APOE {d}") for d in ALL_DIPLOTYPES},
    )


def apoe_genome(rs429358: str | None, rs7412: str | None) -> Genome:
    genome = Genome(source="test")
    if rs429358 is not None:
        genome.add(Genotype(rsid="rs429358", chromosome="19", position=45411941, alleles=rs429358))
    if rs7412 is not None:
        genome.add(Genotype(rsid="rs7412", chromosome="19", position=45412079, alleles=rs7412))
    return genome


# ---------------------------------------------------------------------------
# The diplotype table
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("rs429358", "rs7412", "expected"),
    [
        ("TT", "TT", "e2/e2"),
        ("TT", "CT", "e2/e3"),
        ("TT", "CC", "e3/e3"),  # the common one, ~60% of Europeans
        ("CT", "CC", "e3/e4"),
        ("CC", "CC", "e4/e4"),
        ("CT", "TT", "e1/e2"),
        ("CC", "CT", "e1/e4"),
        ("CC", "TT", "e1/e1"),
    ],
)
def test_unambiguous_diplotypes(rs429358, rs7412, expected):
    finding = derive_apoe(apoe_genome(rs429358, rs7412), apoe_marker())
    assert finding.status is FindingStatus.OK
    assert finding.genotype == expected


def test_allele_order_in_input_does_not_matter():
    assert derive_apoe(apoe_genome("TC", "CC"), apoe_marker()).genotype == "e3/e4"


# ---------------------------------------------------------------------------
# The one genuinely ambiguous call
# ---------------------------------------------------------------------------


def test_double_het_reports_e2_e4_and_says_it_is_ambiguous():
    """rs429358=CT with rs7412=CT is either e2/e4 or e1/e3.

    An unphased array cannot distinguish them. e2/e4 is reported because e1 is
    vanishingly rare, but the ambiguity must reach the user: e2 is protective
    and e4 is the risk allele, so the two readings point opposite ways.
    """
    finding = derive_apoe(apoe_genome("CT", "CT"), apoe_marker())

    assert finding.status is FindingStatus.OK
    assert finding.genotype == "e2/e4"
    assert "e1/e3" in finding.interpretation.detail
    assert "not formally excluded" in finding.interpretation.detail


def test_ambiguity_note_does_not_leak_into_other_diplotypes():
    finding = derive_apoe(apoe_genome("TT", "CC"), apoe_marker())
    assert "e1/e3" not in (finding.interpretation.detail or "")


# ---------------------------------------------------------------------------
# Missing inputs
# ---------------------------------------------------------------------------


def test_missing_rs429358_yields_no_partial_answer():
    """Half an APOE result is worse than none.

    rs7412 alone cannot separate e2 from e1, and 23andMe coverage of rs429358
    has varied by chip version, so this path is routine rather than exotic.
    """
    finding = derive_apoe(apoe_genome(None, "CC"), apoe_marker())
    assert finding.status is FindingStatus.NOT_ON_CHIP
    assert finding.genotype is None


def test_missing_rs7412_yields_no_partial_answer():
    assert derive_apoe(apoe_genome("CT", None), apoe_marker()).status is FindingStatus.NOT_ON_CHIP


def test_no_call_on_either_snp():
    assert derive_apoe(apoe_genome("--", "CC"), apoe_marker()).status is FindingStatus.NO_CALL
    assert derive_apoe(apoe_genome("CT", "--"), apoe_marker()).status is FindingStatus.NO_CALL


def test_impossible_genotype_is_surfaced():
    """APOE SNPs are both C/T. An A or G call means something is wrong upstream."""
    finding = derive_apoe(apoe_genome("AA", "CC"), apoe_marker())
    assert finding.status is FindingStatus.UNKNOWN_GENOTYPE


# ---------------------------------------------------------------------------
# Zygosity and dispatch
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("rs429358", "rs7412", "expected"),
    [
        ("CC", "CC", Zygosity.HOMOZYGOUS_ALT),  # e4/e4
        ("CT", "CC", Zygosity.HETEROZYGOUS),  # e3/e4
        ("TT", "CC", Zygosity.HOMOZYGOUS_REF),  # e3/e3
    ],
)
def test_zygosity_counts_e4_copies(rs429358, rs7412, expected):
    assert derive_apoe(apoe_genome(rs429358, rs7412), apoe_marker()).zygosity is expected


def test_interpret_marker_dispatches_to_derivation():
    """A marker with derivation: apoe must route here, not to genotype lookup."""
    finding = interpret_marker(apoe_genome("TT", "CC"), apoe_marker())
    assert finding.genotype == "e3/e3"


def test_unknown_derivation_name_fails_loudly():
    marker = Marker(
        id="x",
        gene="X",
        name="X",
        category=Category.NEURO,
        evidence=Evidence.WEAK,
        rsids=("rs1", "rs2"),
        derivation="does_not_exist",
        interpretations={"a": Interpretation(summary="a")},
    )
    with pytest.raises(ValueError, match="not registered"):
        interpret_marker(Genome(source="test"), marker)


# ---------------------------------------------------------------------------
# HFE
# ---------------------------------------------------------------------------


def hfe_marker():
    from dnadistiller.catalog import load_catalog

    return next(m for m in load_catalog() if m.id == "hfe_haemochromatosis")


def hfe_genome(c282y: str | None = None, h63d: str | None = None) -> Genome:
    genome = Genome(source="test")
    if c282y is not None:
        genome.add(Genotype(rsid="rs1800562", chromosome="6", position=26092913, alleles=c282y))
    if h63d is not None:
        genome.add(Genotype(rsid="rs1799945", chromosome="6", position=26090951, alleles=h63d))
    return genome


@pytest.mark.parametrize(
    ("c282y", "h63d", "expected"),
    [
        ("GG", "CC", "none"),
        ("GG", "CG", "h63d_carrier"),
        ("GG", "GG", "h63d_homozygous"),
        ("AG", "CC", "c282y_carrier"),
        ("AG", "CG", "compound_heterozygous"),
        ("AA", "CC", "c282y_homozygous"),
        ("AA", "CG", "c282y_homozygous"),
    ],
)
def test_hfe_diplotypes(c282y, h63d, expected):
    from dnadistiller.derive import derive_hfe

    assert derive_hfe(hfe_genome(c282y, h63d), hfe_marker()).genotype == expected


def test_hfe_compound_heterozygote_is_distinguished_from_two_carriers():
    """The reason both SNPs are read together at all.

    One copy of each is a real, milder risk category. Reported as two separate
    markers it would render as two unremarkable carrier results and the category
    would vanish.
    """
    from dnadistiller.derive import derive_hfe

    compound = derive_hfe(hfe_genome("AG", "CG"), hfe_marker())
    carrier = derive_hfe(hfe_genome("AG", "CC"), hfe_marker())

    assert compound.genotype == "compound_heterozygous"
    assert carrier.genotype == "c282y_carrier"
    assert compound.interpretation.summary != carrier.interpretation.summary


def test_hfe_resolves_without_h63d_when_c282y_is_homozygous():
    """A missing H63D is survivable; a missing C282Y is not.

    C282Y homozygosity answers the question on its own, so an absent H63D must
    not throw away a result that matters.
    """
    from dnadistiller.derive import derive_hfe

    finding = derive_hfe(hfe_genome("AA", None), hfe_marker())
    assert finding.status is FindingStatus.OK
    assert finding.genotype == "c282y_homozygous"


def test_hfe_declines_when_h63d_missing_and_c282y_heterozygous():
    """Here H63D decides between carrier and compound het, so it is required."""
    from dnadistiller.derive import derive_hfe

    assert derive_hfe(hfe_genome("AG", None), hfe_marker()).status is FindingStatus.NO_CALL


def test_hfe_declines_without_c282y():
    from dnadistiller.derive import derive_hfe

    assert derive_hfe(hfe_genome(None, "CG"), hfe_marker()).status is FindingStatus.NOT_ON_CHIP


def test_hfe_penetrance_is_stated_not_implied():
    """The wording is the safety feature.

    Most C282Y homozygotes never develop iron overload, and a tool that reports
    the genotype without that number manufactures a patient.
    """
    from dnadistiller.derive import derive_hfe

    detail = derive_hfe(hfe_genome("AA", "CC"), hfe_marker()).interpretation.detail
    assert "28%" in detail
    assert "never become ill" in detail
