"""Tests for selection and redaction.

These are the tests that matter most. Everything else in this project produces a
wrong answer when it breaks; this module produces a privacy breach, and the user
does not find out. So the leak tests below assert on the rendered output rather
than on internal state: what matters is what ends up in the file someone pastes,
not what `redact` believed it was doing.
"""

from __future__ import annotations

import json
import re

import pytest

from dnadistiller.models import (
    Build,
    Category,
    Evidence,
    Finding,
    FindingStatus,
    Genome,
    Interpretation,
    Marker,
    Profile,
    Tier,
    Zygosity,
)
from dnadistiller.output import to_json, to_markdown
from dnadistiller.profile import redact, select_markers

# credentials. The name is the point: if either string reaches output at a tier
# that does not permit it, that is the bug this module exists to catch.
SECRET_RSID = "rs1801133"  # noqa: S105
SECRET_GENOTYPE = "AG"  # noqa: S105


def make_marker(**overrides) -> Marker:
    defaults = dict(
        id="test_marker",
        gene="MTHFR",
        name="Folate enzyme activity",
        category=Category.METHYLATION,
        evidence=Evidence.MODERATE,
        rsids=(SECRET_RSID,),
        interpretations={"AG": Interpretation(summary="Somewhat reduced enzyme activity")},
        effect_allele="A",
        citations=("PMID:12345678",),
        effect_size="OR 1.1",
    )
    defaults.update(overrides)
    return Marker(**defaults)


def make_finding(marker: Marker | None = None) -> Finding:
    marker = marker or make_marker()
    return Finding(
        marker=marker,
        status=FindingStatus.OK,
        genotype=SECRET_GENOTYPE,
        zygosity=Zygosity.HETEROZYGOUS,
        interpretation=marker.interpretations["AG"],
    )


def make_profile(tier: Tier, findings=None) -> Profile:
    return Profile(
        findings=findings if findings is not None else [make_finding()],
        tier=tier,
        source="23andMe",
        build=Build.GRCH37,
    )


# ---------------------------------------------------------------------------
# Leak tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("renderer", [to_markdown, to_json])
def test_minimal_tier_leaks_neither_rsid_nor_genotype_nor_gene(renderer):
    """MINIMAL promises interpretations only. Hold it to that in the actual output."""
    rendered = renderer(make_profile(Tier.MINIMAL))

    assert SECRET_RSID not in rendered
    assert "MTHFR" not in rendered
    # The summary itself must survive, or the tier would be useless.
    assert "Somewhat reduced enzyme activity" in rendered


@pytest.mark.parametrize("renderer", [to_markdown, to_json])
def test_standard_tier_leaks_no_rsid(renderer):
    """STANDARD discloses the gene and copy count, never the rsID or raw call."""
    rendered = renderer(make_profile(Tier.STANDARD))

    assert SECRET_RSID not in rendered
    assert "MTHFR" in rendered
    assert "one copy" in rendered


@pytest.mark.parametrize("renderer", [to_markdown, to_json])
def test_full_tier_discloses_everything(renderer):
    """FULL is the tier that does disclose. Confirm it actually does."""
    rendered = renderer(make_profile(Tier.FULL))

    assert SECRET_RSID in rendered
    assert SECRET_GENOTYPE in rendered
    assert "MTHFR" in rendered


def test_json_never_discloses_more_than_markdown():
    """The renderers must withhold the same things at the same tier.

    The realistic bug is a field added to one renderer and not the other, which
    makes the tier mean whichever format you happened to choose. This caught a
    real one: JSON emitted `trait` at standard tier while Markdown did not.
    """
    for tier in Tier:
        profile = make_profile(tier)
        payload = json.loads(to_json(profile))
        markdown = to_markdown(profile).lower()

        for finding in payload["findings"]:
            for field, value in finding.items():
                if not value:
                    continue
                # `match_method: direct` is intentionally not rendered: it is the
                # unremarkable case and printing it on every line would bury the
                # complement matches that actually warrant a reader's attention.
                if field == "match_method" and value == "direct":
                    continue
                assert str(value).lower() in markdown, (
                    f"{tier}: JSON disclosed {field}={value!r} but Markdown did not. "
                    "The renderers have diverged."
                )


def test_markdown_never_discloses_more_than_json_at_minimal():
    """The reverse direction, which is easy to miss.

    Markdown groups findings under category headings. At minimal tier that would
    disclose the topic of every finding while JSON emits only summaries, so the
    grouping is suppressed. Caught here rather than in review.
    """
    profile = make_profile(Tier.MINIMAL)
    markdown = to_markdown(profile).lower()

    for category in Category:
        assert f"### {category}" not in markdown


def test_redact_refuses_non_reportable_finding():
    """A finding with no interpretation must not render as an empty result."""
    finding = Finding(marker=make_marker(), status=FindingStatus.NOT_ON_CHIP)
    with pytest.raises(ValueError, match="non-reportable"):
        redact(finding, Tier.FULL)


# ---------------------------------------------------------------------------
# Selection
# ---------------------------------------------------------------------------


def test_sensitive_markers_excluded_by_default():
    catalog = [make_marker(), make_marker(id="apoe", gene="APOE", sensitive=True)]
    assert [m.id for m in select_markers(catalog)] == ["test_marker"]


def test_sensitive_markers_included_only_on_request():
    catalog = [make_marker(), make_marker(id="apoe", gene="APOE", sensitive=True)]
    selected = select_markers(catalog, include_sensitive=True)
    assert {m.id for m in selected} == {"test_marker", "apoe"}


def test_weak_markers_excluded_by_default():
    """The default threshold is MODERATE, so weak markers are opt-in."""
    catalog = [make_marker(), make_marker(id="weak_one", evidence=Evidence.WEAK)]
    assert [m.id for m in select_markers(catalog)] == ["test_marker"]


def test_min_evidence_weak_includes_everything():
    catalog = [make_marker(), make_marker(id="weak_one", evidence=Evidence.WEAK)]
    selected = select_markers(catalog, min_evidence=Evidence.WEAK)
    assert len(selected) == 2


def test_min_evidence_strong_excludes_moderate():
    catalog = [make_marker(), make_marker(id="strong_one", evidence=Evidence.STRONG)]
    selected = select_markers(catalog, min_evidence=Evidence.STRONG)
    assert [m.id for m in selected] == ["strong_one"]


def test_category_filter():
    catalog = [make_marker(), make_marker(id="lipid_one", category=Category.LIPIDS)]
    selected = select_markers(catalog, categories={Category.LIPIDS})
    assert [m.id for m in selected] == ["lipid_one"]


def test_sensitive_stays_excluded_even_when_its_category_is_requested():
    """Category and sensitivity are independent gates.

    Asking for --category neuro is not consent to an APOE result; that needs
    --include-sensitive as well. Worth pinning, since making one filter imply
    the other is a tempting convenience.
    """
    catalog = [make_marker(id="apoe", category=Category.NEURO, sensitive=True)]
    assert select_markers(catalog, categories={Category.NEURO}) == []


# ---------------------------------------------------------------------------
# Missing data
# ---------------------------------------------------------------------------


def test_real_catalog_summaries_never_name_their_own_gene():
    """The catalog's own text must respect the tier, not just the renderer.

    `summary` is the only field minimal tier emits. A summary that names its gene
    hands over the gene name no matter how careful the renderer is, and for APOE
    it hands over the genotype too: "one copy of the higher-risk APOE type" IS
    the result. This caught exactly that in the shipped catalog.
    """
    from dnadistiller.catalog import load_catalog

    offenders = [
        (marker.id, genotype)
        for marker in load_catalog()
        for genotype, interp in marker.interpretations.items()
        if marker.gene.upper() in interp.summary.upper()
    ]
    assert not offenders, (
        f"these interpretations name their own gene in `summary`, which minimal "
        f"tier renders verbatim: {offenders}"
    )


def test_real_catalog_summaries_never_contain_an_rsid():
    from dnadistiller.catalog import load_catalog

    offenders = [
        (marker.id, genotype)
        for marker in load_catalog()
        for genotype, interp in marker.interpretations.items()
        if re.search(r"\brs\d+\b", interp.summary)
    ]
    assert not offenders, f"rsID leaked into a summary: {offenders}"


def test_minimal_tier_withholds_gene_names_from_the_not_tested_list():
    """The 'Not tested' section must redact too.

    It is the one place gene names reached the output regardless of tier. A
    promise with an exception for one section is not a promise.
    """
    marker = make_marker(gene="LPA")
    finding = Finding(marker=marker, status=FindingStatus.NOT_ON_CHIP)
    profile = Profile(
        findings=[finding],
        tier=Tier.MINIMAL,
        source="23andMe",
        build=Build.GRCH37,
        missing=[marker],
    )

    rendered = to_markdown(profile)
    assert "LPA" not in rendered
    # The count still has to survive: absence is a result.
    assert "Not tested" in rendered
    assert "(1)" in rendered


def test_untested_markers_are_reported_not_omitted():
    """A silent gap reads as an all-clear. It must be stated instead."""
    marker = make_marker(gene="LPA")
    finding = Finding(marker=marker, status=FindingStatus.NOT_ON_CHIP)
    profile = Profile(
        findings=[finding],
        tier=Tier.STANDARD,
        source="23andMe",
        build=Build.GRCH37,
        missing=[marker],
    )

    rendered = to_markdown(profile)
    assert "Not tested" in rendered
    assert "LPA" in rendered
    assert "not a normal result" in rendered


def test_empty_profile_says_so_rather_than_rendering_blank():
    profile = Profile(findings=[], tier=Tier.STANDARD, source="23andMe", build=Build.GRCH37)
    assert "No reportable findings" in to_markdown(profile)


def test_preamble_is_always_present():
    """The preamble is what keeps an LLM calibrated. It is not optional."""
    for tier in Tier:
        rendered = to_markdown(make_profile(tier))
        assert "not clinically validated" in rendered.lower()
        assert "association is not causation" in rendered.lower()


# ---------------------------------------------------------------------------
# Blood tests
# ---------------------------------------------------------------------------


def test_blood_tests_render_at_every_tier_including_minimal():
    """A test name discloses nothing about you.

    "Ask for a ferritin" is a sentence about medicine, not about this person's
    genome. Withholding it at minimal tier would cost the user the single most
    actionable thing the tool produces and buy no privacy at all.
    """
    from dnadistiller.catalog import load_catalog
    from dnadistiller.interpret import interpret_marker
    from dnadistiller.models import Genotype

    hfe = next(m for m in load_catalog() if m.id == "hfe_haemochromatosis")
    genome = Genome(source="test")
    genome.add(Genotype(rsid="rs1800562", chromosome="6", position=1, alleles="AA"))
    genome.add(Genotype(rsid="rs1799945", chromosome="6", position=2, alleles="CC"))

    for tier in Tier:
        profile = Profile(
            findings=[interpret_marker(genome, hfe)],
            tier=tier,
            source="23andMe",
            build=Build.GRCH37,
        )
        rendered = to_markdown(profile)
        assert "Tests worth asking for" in rendered
        assert "Ferritin" in rendered


def test_blood_test_section_withholds_the_gene_at_minimal():
    """The test name is safe. The gene it came from is not."""
    from dnadistiller.catalog import load_catalog
    from dnadistiller.interpret import interpret_marker
    from dnadistiller.models import Genotype

    hfe = next(m for m in load_catalog() if m.id == "hfe_haemochromatosis")
    genome = Genome(source="test")
    genome.add(Genotype(rsid="rs1800562", chromosome="6", position=1, alleles="AA"))
    genome.add(Genotype(rsid="rs1799945", chromosome="6", position=2, alleles="CC"))

    profile = Profile(
        findings=[interpret_marker(genome, hfe)],
        tier=Tier.MINIMAL,
        source="23andMe",
        build=Build.GRCH37,
    )
    rendered = to_markdown(profile)
    assert "Ferritin" in rendered
    assert "HFE" not in rendered


def test_most_of_the_catalog_admits_a_blood_test_beats_it():
    """The thesis, asserted.

    If this ever flips, the catalog has drifted into pretending genotypes are
    measurements. Genetics is a prior; the assay is the answer.
    """
    from dnadistiller.catalog import load_catalog
    from dnadistiller.models import TestVerdict

    tests = [m.blood_test for m in load_catalog() if m.blood_test]
    supersedes = [t for t in tests if t.verdict is TestVerdict.SUPERSEDES]
    assert len(supersedes) >= len(tests) / 3, (
        "most markers should defer to a measurement; if they do not, check "
        "whether the catalog has started overselling genotypes"
    )
