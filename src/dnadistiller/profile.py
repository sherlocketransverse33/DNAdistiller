"""Selection and redaction: turning findings into something safe to share.

This module is the one that decides what leaves the machine, so it is written
to be read. Two rules shape it:

1. **Opt in, never opt out.** Every filter here starts from nothing and adds
   what was asked for. A bug in an opt-out filter discloses more than the user
   wanted; the same bug in an opt-in filter discloses less. Only one of those
   failure modes is recoverable.

2. **Absence is a result.** A marker the chip never tested is reported as "not
   tested", not omitted. A gap in a list of risks reads as "no risk found",
   which is a claim the data does not support.
"""

from __future__ import annotations

from .models import (
    Category,
    Evidence,
    Finding,
    FindingStatus,
    Genome,
    Marker,
    Profile,
    Severity,
    Tier,
)

#: Ranked weakest-to-strongest so `min_evidence` can be a threshold.
_EVIDENCE_RANK: dict[Evidence, int] = {
    Evidence.WEAK: 0,
    Evidence.MODERATE: 1,
    Evidence.STRONG: 2,
}


def select_markers(
    catalog: list[Marker],
    *,
    categories: set[Category] | None = None,
    min_evidence: Evidence = Evidence.MODERATE,
    include_sensitive: bool = False,
) -> list[Marker]:
    """Choose which markers to even look up.

    Filtering happens here, before interpretation, rather than at render time.
    A marker excluded here is never read out of the genome at all, so it cannot
    reach the output through some later bug — the safe outcome is structural
    rather than something the renderer has to remember to do.

    Args:
        categories: Topics to include. None means all topics — this is the
            catalog-browsing default, not the profile default; `cli.py` passes
            an explicit set.
        min_evidence: Drop markers whose support is weaker than this. Defaults
            to MODERATE, so the weak-but-popular markers (MTHFR being the
            standing example) are opt-in rather than served up alongside real
            findings as though they were equivalent.
        include_sensitive: Whether to include markers flagged `sensitive`,
            notably APOE. Off by default. Genetic counselling treats
            Alzheimer's risk disclosure as something a person decides to
            receive, not something they encounter by accident while reading
            about caffeine.
    """
    threshold = _EVIDENCE_RANK[min_evidence]
    selected = []
    for marker in catalog:
        if categories is not None and marker.category not in categories:
            continue
        if _EVIDENCE_RANK[marker.evidence] < threshold:
            continue
        if marker.sensitive and not include_sensitive:
            continue
        selected.append(marker)
    return selected


def build_profile(
    genome: Genome,
    markers: list[Marker],
    *,
    tier: Tier,
    interpret_fn: object = None,
) -> Profile:
    """Assemble a profile from a genome and a pre-filtered marker list.

    `markers` is expected to have already been through `select_markers`; this
    function does not re-check categories or sensitivity, because having the
    filter in exactly one place is what makes it auditable.

    Markers that could not be read are collected into `Profile.missing` rather
    than dropped, so the renderer can distinguish "we looked and found nothing
    notable" from "we never looked".
    """
    # Imported lazily to keep the module graph acyclic: interpret imports models,
    # and wiring it at the top here would make profile <-> interpret circular.
    if interpret_fn is None:
        from .interpret import interpret_marker

        interpret_fn = interpret_marker

    findings: list[Finding] = []
    missing: list[Marker] = []

    for marker in markers:
        finding = interpret_fn(genome, marker)  # type: ignore[operator]
        if finding.status is FindingStatus.OK:
            findings.append(finding)
        else:
            missing.append(marker)
            findings.append(finding)

    # Blind spots and parse issues carry through to the shared profile. They are
    # facts about the chip and the file, not about the person, so they disclose
    # nothing at any tier — and they are the context most likely to stop a reader
    # treating an untested marker as a negative result.
    return Profile(
        findings=findings,
        tier=tier,
        source=genome.source,
        build=genome.build,
        missing=missing,
        blind_spots=list(genome.blind_spots),
        issues=[i for i in genome.issues if i.severity is not Severity.INFO],
    )


def redact(finding: Finding, tier: Tier) -> dict[str, str]:
    """Render one finding down to the fields its tier permits.

    Returns a plain dict rather than a formatted string so that every output
    format redacts identically — a bug fixed in the Markdown renderer that was
    never fixed in the JSON one is exactly the kind of leak this design is
    meant to make impossible.

    The tiers, concretely::

        FULL      {gene, rsid, genotype, zygosity, summary, evidence, ...}
        STANDARD  {gene, zygosity, summary, evidence}
        MINIMAL   {summary}

    MINIMAL drops the gene name too, which is a deliberate cost. "Reduced
    caffeine metabolism" is discussable; "CYP1A2" is a lookup key into a public
    database, and a list of gene names plus their calls narrows a person down
    almost as efficiently as the rsIDs would.
    """
    if not finding.is_reportable or finding.interpretation is None:
        raise ValueError(f"cannot redact non-reportable finding for {finding.marker.id!r}")

    interp = finding.interpretation
    marker = finding.marker

    if tier is Tier.MINIMAL:
        return {"summary": interp.summary}

    if tier is Tier.STANDARD:
        return {
            "gene": marker.gene,
            "trait": marker.name,
            "zygosity": finding.zygosity.value,
            "summary": interp.summary,
            "evidence": str(marker.evidence),
        }

    return {
        "gene": marker.gene,
        "trait": marker.name,
        "rsid": ", ".join(marker.rsids),
        "genotype": finding.genotype or "",
        "zygosity": finding.zygosity.value,
        "summary": interp.summary,
        "detail": interp.detail,
        "evidence": str(marker.evidence),
        "effect_size": marker.effect_size,
        "category": str(marker.category),
        "citations": "; ".join(marker.citations),
        "ancestry_note": marker.ancestry_note,
        # How the call was matched. Only meaningful at FULL, where the rsID and
        # genotype are present and a reader could check them against dbSNP; the
        # answer to "why does your G/A file show as C/T" lives here.
        "match_method": str(finding.match_method),
    }
