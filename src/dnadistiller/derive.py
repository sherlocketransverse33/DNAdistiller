"""Markers that need logic rather than a genotype lookup.

Most markers are one SNP and a dictionary. A few are not: APOE's clinically
meaningful unit is an ε2/ε3/ε4 haplotype spread across two SNPs, and neither
SNP means anything on its own. Handlers live here and are named from a marker's
`derivation:` field in YAML.

Keep this module small. Every entry is bespoke logic that a catalog reviewer
cannot check by reading YAML, so a marker belongs here only when it genuinely
cannot be expressed as data.
"""

from __future__ import annotations

from collections.abc import Callable

from .models import (
    Finding,
    FindingStatus,
    Genome,
    Marker,
    MatchMethod,
    Zygosity,
)

# ---------------------------------------------------------------------------
# APOE
# ---------------------------------------------------------------------------

#: The two SNPs whose combination defines the APOE isoform.
_APOE_RS429358 = "rs429358"
_APOE_RS7412 = "rs7412"

#: Haplotype definitions, on the GRCh37 plus strand — the strand 23andMe and
#: AncestryDNA both report on:
#:
#:     ε2 = rs429358(T) + rs7412(T)
#:     ε3 = rs429358(T) + rs7412(C)     <- the common one
#:     ε4 = rs429358(C) + rs7412(C)
#:     ε1 = rs429358(C) + rs7412(T)     <- vanishingly rare
#:
#: An array reports unphased genotypes: it says "one C and one T at this
#: position" without saying which copy of chromosome 19 each sits on. For eight
#: of the nine genotype pairs only one diplotype is consistent, so phase does
#: not matter. The ninth is genuinely ambiguous — see below.
_APOE_DIPLOTYPES: dict[tuple[str, str], str] = {
    ("TT", "TT"): "e2/e2",
    ("TT", "CT"): "e2/e3",
    ("TT", "CC"): "e3/e3",
    ("CT", "CC"): "e3/e4",
    ("CC", "CC"): "e4/e4",
    # ε1-containing diplotypes. Real, but rare enough that most people will
    # never see one; resolved rather than dropped so the output says something
    # true instead of "unrecognised".
    ("CT", "TT"): "e1/e2",
    ("CC", "CT"): "e1/e4",
    ("CC", "TT"): "e1/e1",
}

#: The one pair phase cannot resolve. rs429358=CT with rs7412=CT is consistent
#: with either ε2/ε4 or ε1/ε3, and an unphased array cannot distinguish them.
#: ε1 is so rare — a handful of families in the literature — that ε2/ε4 is the
#: near-certain answer, and reporting it is standard practice. The ambiguity is
#: still surfaced in `detail` rather than hidden, because the two have opposite
#: implications: ε2 is protective and ε4 is the risk allele, so a person told
#: "ε2/ε4" deserves to know the call rests on a frequency argument.
_APOE_AMBIGUOUS = ("CT", "CT")
_APOE_AMBIGUOUS_RESULT = "e2/e4"
_APOE_AMBIGUOUS_NOTE = (
    "This genotype is consistent with either e2/e4 or e1/e3. An array cannot tell "
    "them apart because it does not resolve which chromosome each allele sits on. "
    "e2/e4 is reported because e1 is extraordinarily rare; e1/e3 is not formally "
    "excluded. Confirm with a clinical-grade test before acting on this."
)


def derive_apoe(genome: Genome, marker: Marker) -> Finding:
    """Resolve APOE ε-genotype from rs429358 and rs7412.

    Requires both SNPs. If either is missing the result is NOT_ON_CHIP rather
    than a partial answer — a lone rs7412 call cannot distinguish ε2 from ε1,
    and reporting half an APOE result would be worse than reporting none.

    23andMe's coverage of rs429358 has varied across chip versions, so a missing
    result here is a routine outcome and not a failure.
    """
    first = genome.get(_APOE_RS429358)
    second = genome.get(_APOE_RS7412)

    if first is None or second is None:
        return Finding(marker=marker, status=FindingStatus.NOT_ON_CHIP)

    if first.is_no_call or second.is_no_call:
        return Finding(marker=marker, status=FindingStatus.NO_CALL)

    key = (first.canonical(), second.canonical())

    # Both SNPs are C/T sites, so the complement of a valid genotype is another
    # valid genotype and a flip cannot be detected from the call alone. Both
    # major providers report APOE on the plus strand, so we match as-read and
    # decline rather than guess if that fails.
    if key == _APOE_AMBIGUOUS:
        diplotype = _APOE_AMBIGUOUS_RESULT
        ambiguous = True
    elif key in _APOE_DIPLOTYPES:
        diplotype = _APOE_DIPLOTYPES[key]
        ambiguous = False
    else:
        return Finding(
            marker=marker,
            status=FindingStatus.UNKNOWN_GENOTYPE,
            genotype=f"{_APOE_RS429358}={key[0]}, {_APOE_RS7412}={key[1]}",
        )

    interpretation = marker.interpretations.get(diplotype)
    if interpretation is None:
        return Finding(marker=marker, status=FindingStatus.UNKNOWN_GENOTYPE, genotype=diplotype)

    if ambiguous:
        detail = f"{interpretation.detail} {_APOE_AMBIGUOUS_NOTE}".strip()
        interpretation = type(interpretation)(
            summary=interpretation.summary,
            detail=detail,
            population_frequency=interpretation.population_frequency,
        )

    return Finding(
        marker=marker,
        status=FindingStatus.OK,
        genotype=diplotype,
        zygosity=_apoe_zygosity(diplotype),
        interpretation=interpretation,
        match_method=MatchMethod.DERIVED,
    )


def _apoe_zygosity(diplotype: str) -> Zygosity:
    """Copies of ε4, which is the axis people mean when they ask about APOE."""
    return {
        0: Zygosity.HOMOZYGOUS_REF,
        1: Zygosity.HETEROZYGOUS,
        2: Zygosity.HOMOZYGOUS_ALT,
    }[diplotype.count("e4")]


# ---------------------------------------------------------------------------
# HFE
# ---------------------------------------------------------------------------

_HFE_C282Y = "rs1800562"
_HFE_H63D = "rs1799945"


def derive_hfe(genome: Genome, marker: Marker) -> Finding:
    """Resolve HFE genotype from C282Y and H63D together.

    Neither SNP is worth reporting alone. C282Y carries essentially all the risk;
    H63D on its own is close to null and reporting it as a finding would alarm a
    quarter of Europeans over nothing. The reason to read them together is the
    compound heterozygote, one copy of each, which is a real if modest risk
    category that two separate results would never reveal.

    Phase would normally be a problem here, since an array does not say which
    chromosome each allele sits on, and "one of each" only means anything if they
    are on opposite chromosomes. It is safe in this specific case: C282Y and H63D
    are in strong negative linkage disequilibrium and essentially never occur
    together in cis, so a double heterozygote is a compound heterozygote in
    trans with high confidence. That is a fact about this locus, not a general
    licence to infer phase.

    Unlike APOE, a missing H63D is survivable: C282Y alone still answers the
    question that matters. Only a missing C282Y makes the result meaningless.
    """
    c282y = genome.get(_HFE_C282Y)
    h63d = genome.get(_HFE_H63D)

    if c282y is None:
        return Finding(marker=marker, status=FindingStatus.NOT_ON_CHIP)
    if c282y.is_no_call:
        return Finding(marker=marker, status=FindingStatus.NO_CALL)

    # Plus strand: rs1800562 G>A with A the risk allele, rs1799945 C>G with G
    # the minor allele. Both sites are unambiguous (G/A and C/G are... note C/G
    # IS palindromic, so H63D must never be complement-matched; we read it as
    # written, which is correct for every provider that reports plus strand).
    c282y_count = c282y.canonical().count("A")

    if c282y_count == 2:
        return _hfe_finding(marker, "c282y_homozygous")

    # H63D only changes the answer for C282Y heterozygotes, so a missing or
    # failed H63D is only fatal to the result in that one case.
    h63d_count = 0
    if h63d is not None and not h63d.is_no_call:
        h63d_count = h63d.canonical().count("G")
    elif c282y_count == 1:
        return Finding(marker=marker, status=FindingStatus.NO_CALL)

    if c282y_count == 1:
        return _hfe_finding(marker, "compound_heterozygous" if h63d_count >= 1 else "c282y_carrier")

    if h63d_count == 2:
        return _hfe_finding(marker, "h63d_homozygous")
    if h63d_count == 1:
        return _hfe_finding(marker, "h63d_carrier")
    return _hfe_finding(marker, "none")


def _hfe_finding(marker: Marker, key: str) -> Finding:
    interpretation = marker.interpretations.get(key)
    if interpretation is None:
        return Finding(marker=marker, status=FindingStatus.UNKNOWN_GENOTYPE, genotype=key)

    # Zygosity here counts copies of C282Y, the allele that carries the risk.
    zygosity = {
        "c282y_homozygous": Zygosity.HOMOZYGOUS_ALT,
        "c282y_carrier": Zygosity.HETEROZYGOUS,
        "compound_heterozygous": Zygosity.HETEROZYGOUS,
    }.get(key, Zygosity.HOMOZYGOUS_REF)

    return Finding(
        marker=marker,
        status=FindingStatus.OK,
        genotype=key,
        zygosity=zygosity,
        interpretation=interpretation,
        match_method=MatchMethod.DERIVED,
    )


#: Registry consulted by `interpret.interpret_marker` via a marker's
#: `derivation:` field.
DERIVATIONS: dict[str, Callable[[Genome, Marker], Finding]] = {
    "apoe": derive_apoe,
    "hfe": derive_hfe,
}
