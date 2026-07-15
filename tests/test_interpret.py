"""Tests for reading a genome against the catalog.

Concentrated on strand handling, which is the most common way a tool in this
space reports a confidently wrong genotype, and on the difference between the
several ways a marker can fail to produce a result.
"""

from __future__ import annotations

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


def marker(**overrides) -> Marker:
    defaults = dict(
        id="m",
        gene="GENE",
        name="Trait",
        category=Category.METABOLIC,
        evidence=Evidence.MODERATE,
        rsids=("rs123",),
        citations=("PMID:1",),
        interpretations={
            "GG": Interpretation(summary="two G"),
            "AG": Interpretation(summary="one A"),
            "AA": Interpretation(summary="two A"),
        },
        effect_allele="A",
    )
    defaults.update(overrides)
    return Marker(**defaults)


def genome_with(alleles: str, rsid: str = "rs123", chromosome: str = "1") -> Genome:
    g = Genome(source="test")
    g.add(Genotype(rsid=rsid, chromosome=chromosome, position=1, alleles=alleles))
    return g


# ---------------------------------------------------------------------------
# Straightforward matching
# ---------------------------------------------------------------------------


def test_direct_match():
    finding = interpret_marker(genome_with("AG"), marker())
    assert finding.status is FindingStatus.OK
    assert finding.interpretation.summary == "one A"
    assert finding.strand_flipped is False


def test_allele_order_does_not_matter():
    """The array reports an unordered pair, so GA and AG are the same call."""
    assert interpret_marker(genome_with("GA"), marker()).interpretation.summary == "one A"


def test_lowercase_rsid_lookup():
    genome = genome_with("AG", rsid="RS123")
    assert interpret_marker(genome, marker()).status is FindingStatus.OK


# ---------------------------------------------------------------------------
# Strand handling
# ---------------------------------------------------------------------------


def test_strand_flip_is_applied_when_direct_match_fails():
    """A catalog on the minus strand still resolves.

    Catalog is keyed A/G. A file reporting the same call on the other strand
    says T/C. That must resolve to the same interpretation, flagged as flipped.
    """
    finding = interpret_marker(genome_with("CT"), marker())
    assert finding.status is FindingStatus.OK
    assert finding.interpretation.summary == "one A"
    assert finding.strand_flipped is True


def test_direct_match_wins_over_flip():
    """A genotype valid on both strands must not be flipped.

    Catalog has GG and AA. A CC call complements to GG, but AA is a direct hit
    and must be preferred. Flipping first would corrupt every correct marker.
    """
    finding = interpret_marker(genome_with("AA"), marker())
    assert finding.interpretation.summary == "two A"
    assert finding.strand_flipped is False


def test_strand_ambiguous_marker_refuses_rather_than_guesses():
    """At an A/T site the complement of a valid genotype is also valid.

    There is no way to tell a flip from a real call, so guessing has a 50%
    chance of inverting the result. Reporting it unread is the honest outcome.
    """
    ambiguous = marker(
        interpretations={
            "AA": Interpretation(summary="two A"),
            "AT": Interpretation(summary="one A"),
            "TT": Interpretation(summary="two T"),
        },
        effect_allele="A",
    )
    finding = interpret_marker(genome_with("CG"), ambiguous)
    assert finding.status is FindingStatus.UNKNOWN_GENOTYPE


def test_cg_site_is_also_treated_as_ambiguous():
    ambiguous = marker(
        interpretations={
            "CC": Interpretation(summary="two C"),
            "CG": Interpretation(summary="one C"),
            "GG": Interpretation(summary="two G"),
        },
        effect_allele="C",
    )
    finding = interpret_marker(genome_with("AT"), ambiguous)
    assert finding.status is FindingStatus.UNKNOWN_GENOTYPE


# ---------------------------------------------------------------------------
# Failure modes, which are distinct and must stay distinct
# ---------------------------------------------------------------------------


def test_marker_absent_from_chip():
    assert interpret_marker(Genome(source="test"), marker()).status is FindingStatus.NOT_ON_CHIP


def test_no_call_is_distinct_from_not_on_chip():
    """Tested-but-failed and never-tested are different facts about the data."""
    assert interpret_marker(genome_with("--"), marker()).status is FindingStatus.NO_CALL


def test_unrecognised_genotype_is_surfaced_not_dropped():
    """A call valid on neither strand is our bug or a build mismatch, not a result.

    AC is the test case rather than the more obvious CC: at an A/G site, CC is
    the legitimate minus-strand reading of GG and correctly resolves. AC has no
    reading on either strand (it complements to GT), so it must be surfaced.
    """
    finding = interpret_marker(genome_with("AC"), marker())
    assert finding.status is FindingStatus.UNKNOWN_GENOTYPE
    assert finding.genotype == "AC"


def test_minus_strand_homozygote_resolves_rather_than_erroring():
    """At an A/G site, CC is minus-strand GG. It is data, not a third allele."""
    finding = interpret_marker(genome_with("CC"), marker())
    assert finding.status is FindingStatus.OK
    assert finding.interpretation.summary == "two G"
    assert finding.strand_flipped is True


def test_no_call_never_produces_an_interpretation():
    assert interpret_marker(genome_with("--"), marker()).interpretation is None


# ---------------------------------------------------------------------------
# Zygosity
# ---------------------------------------------------------------------------


def test_zygosity_counts_effect_allele():
    assert interpret_marker(genome_with("AA"), marker()).zygosity is Zygosity.HOMOZYGOUS_ALT
    assert interpret_marker(genome_with("AG"), marker()).zygosity is Zygosity.HETEROZYGOUS
    assert interpret_marker(genome_with("GG"), marker()).zygosity is Zygosity.HOMOZYGOUS_REF


def test_zygosity_counted_after_strand_flip():
    """A flipped call must be counted on the catalog's strand, not the file's.

    CT flips to AG: one copy of the A effect allele. Counting before the flip
    would find zero, and report the opposite of the truth.
    """
    assert interpret_marker(genome_with("CT"), marker()).zygosity is Zygosity.HETEROZYGOUS


def test_zygosity_unknown_when_marker_declares_no_effect_allele():
    assert interpret_marker(genome_with("AG"), marker(effect_allele=None)).zygosity is (
        Zygosity.UNKNOWN
    )


# ---------------------------------------------------------------------------
# Match provenance
# ---------------------------------------------------------------------------


def test_direct_match_is_recorded_as_direct():
    from dnadistiller.models import MatchMethod

    assert interpret_marker(genome_with("AG"), marker()).match_method is MatchMethod.DIRECT


def test_complement_match_is_recorded_not_hidden():
    """A complement match is an inference, and it must say so.

    The file said CT; we report AG. Anyone checking that genotype against dbSNP
    sees letters that do not agree with their own file and concludes the tool is
    broken, unless the profile tells them a strand flip happened. Recording the
    method is what makes the reasoning inspectable instead of magic.
    """
    from dnadistiller.models import MatchMethod

    finding = interpret_marker(genome_with("CT"), marker())
    assert finding.match_method is MatchMethod.COMPLEMENT
    assert finding.strand_flipped is True
