"""Tests for the catalog validator.

The centrepiece is `test_the_tpmt_inversion_is_rejected`. It reconstructs a real
defect that shipped in another consumer-DNA tool and asserts our validator would
have stopped it. That bug is the best argument this project has for why marker
data needs a schema and a build-time check rather than careful reviewers.
"""

from __future__ import annotations

import pytest

from dnadistiller.catalog.validate import (
    check_citations,
    check_frequencies_are_plausible,
    check_interpretation_coverage,
    check_marker_is_not_inverted,
    check_strong_claims,
)
from dnadistiller.models import Category, Evidence, Interpretation, Marker


def marker(**overrides) -> Marker:
    defaults = dict(
        id="m",
        gene="GENE",
        name="Trait",
        category=Category.PHARMACOGENOMICS,
        evidence=Evidence.MODERATE,
        rsids=("rs1",),
        citations=("PMID:1",),
        ancestry_note="European cohorts.",
        interpretations={
            "TT": Interpretation(summary="Typical.", population_frequency=0.92),
            "CT": Interpretation(summary="Intermediate.", population_frequency=0.08),
            "CC": Interpretation(summary="Reduced.", population_frequency=0.002),
        },
        effect_allele="C",
    )
    defaults.update(overrides)
    return Marker(**defaults)


# ---------------------------------------------------------------------------
# The inversion check
# ---------------------------------------------------------------------------


def test_correctly_oriented_marker_passes():
    assert check_marker_is_not_inverted(marker()) == []


def test_the_tpmt_inversion_is_rejected():
    """Reconstructs a real defect that shipped, published, and passed its own tests.

    Another tool's catalog declared TPMT rs1142345 as::

        "TT": risk    -> "Poor metabolizer. Standard thiopurine doses can cause
                          fatal myelosuppression."
        "CC": normal  -> "Normal metabolizer."

    `T` is the reference allele at rs1142345, carried by about 96% of people; `C`
    at about 4% is the *3C variant. The entry is exactly backwards, so nearly
    every user who ran it was told they risked a fatal drug reaction. It reached
    their published sample report, and their test suite asserted the wrong tier,
    which is what allowed it to survive review.

    Reviewers cannot reliably catch this. The entry looks completely plausible;
    you have to already know which allele is the reference. A machine comparing
    the declared frequency against the declared effect allele catches it every
    time, which is why this check exists.
    """
    inverted = marker(
        gene="TPMT",
        rsids=("rs1142345",),
        # Effect allele declared as the reference allele: the inversion.
        effect_allele="T",
        interpretations={
            "TT": Interpretation(
                summary="Poor metaboliser. Standard thiopurine doses can cause "
                "fatal myelosuppression.",
                population_frequency=0.957,
            ),
            "CT": Interpretation(summary="Intermediate metaboliser.", population_frequency=0.04),
            "CC": Interpretation(summary="Normal metaboliser.", population_frequency=0.002),
        },
    )

    problems = check_marker_is_not_inverted(inverted)
    assert problems, "the validator failed to catch an inverted marker"
    assert "inverted" in problems[0]
    assert "96%" in problems[0]


def test_inversion_check_is_silent_without_a_declared_frequency():
    """No frequency, no opinion. The check must not invent one."""
    quiet = marker(
        interpretations={
            "TT": Interpretation(summary="Typical."),
            "CT": Interpretation(summary="Intermediate."),
            "CC": Interpretation(summary="Reduced."),
        }
    )
    assert check_marker_is_not_inverted(quiet) == []


def test_inversion_check_is_silent_without_an_effect_allele():
    assert check_marker_is_not_inverted(marker(effect_allele=None)) == []


# ---------------------------------------------------------------------------
# Frequency plausibility
# ---------------------------------------------------------------------------


def test_plausible_frequencies_pass():
    assert check_frequencies_are_plausible(marker()) == []


def test_percentages_written_as_proportions_are_caught():
    """The order-of-magnitude slip: 45 where 0.45 was meant."""
    wrong = marker(
        interpretations={
            "TT": Interpretation(summary="a", population_frequency=0.92),
            "CT": Interpretation(summary="b", population_frequency=0.9),
            "CC": Interpretation(summary="c", population_frequency=0.9),
        }
    )
    problems = check_frequencies_are_plausible(wrong)
    assert problems
    assert "sum to 2.72" in problems[0]


def test_frequencies_that_undershoot_are_caught():
    sparse = marker(
        interpretations={
            "TT": Interpretation(summary="a", population_frequency=0.1),
            "CT": Interpretation(summary="b", population_frequency=0.05),
            "CC": Interpretation(summary="c", population_frequency=0.01),
        }
    )
    assert check_frequencies_are_plausible(sparse)


# ---------------------------------------------------------------------------
# The other guards
# ---------------------------------------------------------------------------


def test_missing_common_genotype_is_caught():
    """A catalog describing only the risk genotype makes every profile bad news."""
    partial = marker(
        interpretations={
            "CT": Interpretation(summary="Intermediate."),
            "CC": Interpretation(summary="Reduced."),
        }
    )
    problems = check_interpretation_coverage(partial)
    assert problems
    assert "TT" in problems[0]


@pytest.mark.parametrize("citation", ["Smith et al 2019", "https://example.com/paper", "12345678"])
def test_unresolvable_citations_are_rejected(citation):
    problems = check_citations(marker(citations=(citation,)))
    assert problems
    assert "PMID" in problems[0]


@pytest.mark.parametrize("citation", ["PMID:16522833", "DOI:10.1073/pnas.0801030105"])
def test_resolvable_citations_pass(citation):
    assert check_citations(marker(citations=(citation,))) == []


def test_strong_claims_must_state_an_effect_size():
    problems = check_strong_claims(marker(evidence=Evidence.STRONG, effect_size=""))
    assert any("effect_size" in p for p in problems)


def test_strong_claims_must_state_ancestry():
    problems = check_strong_claims(
        marker(evidence=Evidence.STRONG, effect_size="OR 1.4", ancestry_note="")
    )
    assert any("ancestry_note" in p for p in problems)


# ---------------------------------------------------------------------------
# The shipped catalog
# ---------------------------------------------------------------------------


def test_the_shipped_catalog_passes_every_check():
    """What CI runs, asserted here so a bad marker fails the test suite too."""
    from dnadistiller.catalog import load_catalog
    from dnadistiller.catalog.validate import PER_MARKER_CHECKS, check_duplicate_rsids

    catalog = load_catalog()
    problems = [
        f"{m.id}: {problem}" for m in catalog for check in PER_MARKER_CHECKS for problem in check(m)
    ]
    problems.extend(check_duplicate_rsids(catalog))
    assert not problems, problems
