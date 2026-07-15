"""Rendering a profile for sharing.

Every renderer goes through `profile.redact` rather than reading `Finding`
directly. That is the whole discipline of this package: the redaction rules live
in one function, and a renderer physically cannot disclose a field that
`redact` did not hand it. A Markdown renderer that respects the tier and a JSON
renderer that quietly dumps the full object is the obvious way to build this
tool wrong.
"""

from __future__ import annotations

import json
from collections import defaultdict

from ..models import (
    BloodTest,
    Category,
    Finding,
    FindingStatus,
    Marker,
    Profile,
    TestVerdict,
    Tier,
)
from ..profile import redact
from ..prompts import build_preamble, suggest_prompts

_HEADER_NOTE = (
    "Generated locally by DNAdistiller (https://github.com/Lucaks3/DNAdistiller). "
    "No data left this machine during generation."
)


def to_markdown(profile: Profile, *, include_prompts: bool = True) -> str:
    """Render a profile as Markdown, ready to paste into a chat window.

    Markdown is the default because it survives a paste into an LLM intact and
    stays readable to the person deciding whether to send it — which is the
    review step the whole tiering system depends on. A format you cannot skim
    is one you cannot meaningfully consent to sharing.
    """
    lines: list[str] = ["# Longevity genetic profile", ""]
    lines.append(build_preamble(profile))
    lines.append("")

    reportable = profile.reportable
    if not reportable:
        lines.append(
            "No reportable findings. Every requested marker was either absent from "
            "this chip or failed to produce a call. See 'Not tested' below."
        )
    elif profile.tier is Tier.MINIMAL:
        # Flat list, no category headings. Grouping under "Methylation" would
        # disclose the topic of each finding, which MINIMAL does not promise and
        # the JSON renderer does not emit. The two formats have to withhold the
        # same things or the tier means whichever one you happened to pick.
        lines.append("## Findings")
        lines.append("")
        for finding in sorted(reportable, key=_summary_of):
            lines.extend(_render_finding(finding, profile.tier))
        lines.append("")
    else:
        by_category: dict[Category, list[Finding]] = defaultdict(list)
        for finding in reportable:
            by_category[finding.marker.category].append(finding)

        lines.append("## Findings")
        lines.append("")
        for category in sorted(by_category, key=str):
            lines.append(f"### {str(category).replace('_', ' ').title()}")
            lines.append("")
            for finding in sorted(by_category[category], key=lambda f: f.marker.gene):
                lines.extend(_render_finding(finding, profile.tier))
            lines.append("")

    lines.extend(_render_blood_tests(profile))
    lines.extend(_render_missing(profile))
    lines.extend(_render_blind_spots(profile))
    lines.extend(_render_provenance(profile))

    if include_prompts:
        lines.extend(_render_prompts(profile))

    return "\n".join(lines).rstrip() + "\n"


def _summary_of(finding: Finding) -> str:
    """Sort key for MINIMAL tier, where the summary is all there is to sort on."""
    return finding.interpretation.summary if finding.interpretation else ""


def _render_finding(finding: Finding, tier: Tier) -> list[str]:
    fields = redact(finding, tier)

    if tier is Tier.MINIMAL:
        return [f"- {fields['summary']}"]

    if tier is Tier.STANDARD:
        line = f"- **{fields['gene']}** ({fields['trait']}): {fields['summary']}"
        line += f" ({fields['zygosity']}, evidence: {fields['evidence']})"
        return [line]

    lines = [
        f"- **{fields['gene']}** ({fields['rsid']}): `{fields['genotype']}` ({fields['zygosity']})"
    ]
    lines.append(f"  - Trait: {fields['trait']}")
    lines.append(f"  - {fields['summary']}")
    if fields["detail"]:
        lines.append(f"  - {fields['detail']}")
    if fields["effect_size"]:
        lines.append(f"  - Effect size: {fields['effect_size']}")
    lines.append(f"  - Evidence: {fields['evidence']}")
    if fields["ancestry_note"]:
        lines.append(f"  - Ancestry: {fields['ancestry_note']}")
    if fields["citations"]:
        lines.append(f"  - Citations: {fields['citations']}")

    # Only worth saying when it was not a plain direct match. A complement match
    # means we inferred the file reported on the other strand and rewrote the
    # call to suit; anyone checking this genotype against dbSNP needs to know
    # that, or the letters will not agree and they will assume we are wrong.
    if fields["match_method"] != "direct":
        lines.append(f"  - Matched by: {fields['match_method']}")
    return lines


def _render_missing(profile: Profile) -> list[str]:
    """List markers that could not be read.

    Reported rather than omitted. A list of risk findings with silent gaps reads
    as an all-clear on everything absent, which is a claim an array cannot
    support — it tests preselected positions and says nothing about the rest.

    Gene names are withheld at MINIMAL tier, which costs some usefulness and is
    still correct: the tier promises no gene names, and a promise with an
    exception for one section is not a promise. Counts alone still carry the
    point, which is that the gaps are gaps rather than clean results.
    """
    if not profile.missing:
        return []

    not_on_chip = [f for f in profile.findings if f.status is FindingStatus.NOT_ON_CHIP]
    no_call = [f for f in profile.findings if f.status is FindingStatus.NO_CALL]
    unknown = [f for f in profile.findings if f.status is FindingStatus.UNKNOWN_GENOTYPE]

    named = profile.tier is not Tier.MINIMAL

    def label(findings: list[Finding]) -> str:
        if not named:
            return ""
        return ": " + ", ".join(sorted(f.marker.gene for f in findings))

    lines = ["## Not tested", "", "These were requested but could not be read. "]
    lines.append("Absence here means no result, not a normal result.")
    lines.append("")

    if not_on_chip:
        lines.append(f"- **Not on this chip** ({len(not_on_chip)}){label(not_on_chip)}")
        lines.append(
            "  - This provider's array never tested these positions. A different "
            "provider or chip version may include them."
        )
    if no_call:
        lines.append(f"- **No call** ({len(no_call)}){label(no_call)}")
        lines.append("  - Tested, but the assay did not return a result.")
    if unknown:
        lines.append(f"- **Unrecognised genotype** ({len(unknown)}){label(unknown)}")
        lines.append(
            "  - A call came back that the catalog does not know how to read. This is "
            "more likely our bug than your biology. Please report it."
        )

    lines.append("")
    return lines


def _render_blind_spots(profile: Profile) -> list[str]:
    """State what the source array could not see, at every tier.

    Aimed at whoever reads the profile next, usually a language model. Handed a
    list of genetic findings with no framing, a model treats the absences as
    negatives and the small effects as predictions. This section is the cheapest
    available correction to the first half of that, and it costs no privacy: it
    describes the chip, not the person.
    """
    if not profile.blind_spots and not profile.issues:
        return []

    lines = ["## What this file cannot show", ""]

    for issue in profile.issues:
        lines.append(f"- **{issue.severity.value.upper()}:** {issue.message}")
    for spot in profile.blind_spots:
        lines.append(f"- {spot}")

    lines.append("")
    return lines


def _render_blood_tests(profile: Profile) -> list[str]:
    """The measurements worth having, derived from what was actually found.

    Rendered at every tier, including MINIMAL, because a test name discloses
    nothing about you. "Ask for a ferritin test" is a sentence about medicine,
    not about your genome, and anyone could have said it.

    Split by whether the test supersedes the genotype, because that distinction
    is the honest core of the whole tool. For most of this catalog the blood test
    is simply better: it measures what the variant only predicts, it costs about
    as much as lunch, and it is about you today rather than about a population.
    Saying so costs us the mystique and buys the credibility.
    """
    tests: list[tuple[Marker, BloodTest]] = [
        (f.marker, f.marker.blood_test) for f in profile.reportable if f.marker.blood_test
    ]
    if not tests:
        return []

    def of(verdict: TestVerdict) -> list[tuple[Marker, BloodTest]]:
        return [(m, t) for m, t in tests if t.verdict is verdict]

    lines = ["## Tests worth asking for", ""]
    lines.append(
        "Genetics is a prior. These are the measurements. Most of them cost very "
        "little and tell you about your body now rather than about a population."
    )
    lines.append("")

    named = profile.tier is not Tier.MINIMAL

    def render_named(pairs: list[tuple[Marker, BloodTest]]) -> None:
        for marker, test in sorted(pairs, key=lambda pair: pair[1].name):
            label = f" ({marker.gene})" if named else ""
            lines.append(f"- **{test.name}**{label}")
            if test.note:
                lines.append(f"  - {test.note}")

    if better := of(TestVerdict.SUPERSEDES):
        lines.append("**These beat the genotype. Get the number instead.**")
        lines.append("")
        render_named(better)
        lines.append("")

    if also := of(TestVerdict.COMPLEMENTS):
        lines.append("**Worth having alongside:**")
        lines.append("")
        render_named(also)
        lines.append("")

    if no_test := of(TestVerdict.NONE_EXISTS):
        # The heading carries the point: these are the markers where knowing the
        # genotype is the whole value, because nothing measures them in time.
        lines.append("**No test warns you in time. For these, the genotype is the point.**")
        lines.append("")
        for marker, test in sorted(no_test, key=lambda pair: pair[0].gene):
            label = f"**{marker.gene}**" if named else "**One finding**"
            lines.append(f"- {label}: {test.note}")
        lines.append("")

    return lines


def _render_provenance(profile: Profile) -> list[str]:
    return [
        "## About this profile",
        "",
        f"- Source: {profile.source} export"
        + (f", reference build {profile.build.value}" if profile.build.value else ""),
        f"- Disclosure tier: {profile.tier}",
        f"- Markers reported: {len(profile.reportable)} of {len(profile.findings)} requested",
        f"- {_HEADER_NOTE}",
        "",
    ]


def _render_prompts(profile: Profile) -> list[str]:
    prompts = suggest_prompts(profile)
    if not prompts:
        return []

    lines = [
        "---",
        "",
        "## Questions worth asking",
        "",
        "Suggested follow-ups for whoever is reading this with you:",
        "",
    ]
    lines.extend(f"{i}. {prompt}" for i, prompt in enumerate(prompts, start=1))
    lines.append("")
    return lines


def to_json(profile: Profile, *, include_prompts: bool = True) -> str:
    """Render a profile as JSON, for piping into other tools.

    Redacts through the same path as Markdown — see the module docstring.
    """
    payload: dict[str, object] = {
        "schema": "dnadistiller/profile/v1",
        "tier": str(profile.tier),
        "source": profile.source,
        "build": profile.build.value or None,
        "preamble": build_preamble(profile),
        "findings": [redact(f, profile.tier) for f in profile.reportable],
        "recommended_tests": [
            {
                "test": f.marker.blood_test.name,
                "note": f.marker.blood_test.note,
                "verdict": str(f.marker.blood_test.verdict),
                **({"gene": f.marker.gene} if profile.tier is not Tier.MINIMAL else {}),
            }
            for f in profile.reportable
            if f.marker.blood_test
        ],
        "blind_spots": profile.blind_spots,
        "parse_warnings": [
            {"severity": i.severity.value, "code": i.code, "message": i.message}
            for i in profile.issues
        ],
        "not_tested": [
            {
                "gene": f.marker.gene,
                "trait": f.marker.name,
                "reason": f.status.value,
            }
            for f in profile.findings
            if f.status is not FindingStatus.OK
        ],
    }

    if include_prompts:
        payload["suggested_prompts"] = suggest_prompts(profile)

    return json.dumps(payload, indent=2, ensure_ascii=False) + "\n"


RENDERERS = {
    "md": to_markdown,
    "markdown": to_markdown,
    "json": to_json,
}

__all__ = ["RENDERERS", "to_json", "to_markdown"]
