"""Reading a genome against the marker catalog.

The interesting problem here is not lookup, it is knowing when *not* to answer.
An array can fail to call a position, a chip can omit it entirely, and a call
can come back in a form the catalog does not recognise. Each of those is a
different thing to tell the user, and none of them is "no risk found".

The other trap is strand orientation, handled in `_match`. It is the most
common way tools in this space produce confidently wrong output.
"""

from __future__ import annotations

from .models import (
    Finding,
    FindingStatus,
    Genome,
    Genotype,
    Interpretation,
    Marker,
    MatchMethod,
    Zygosity,
)


def interpret_marker(genome: Genome, marker: Marker) -> Finding:
    """Read one marker out of a genome.

    Never raises for missing or malformed data: a marker that cannot be read
    yields a Finding carrying the reason. Callers render that reason; an
    exception here would tempt a caller into a bare `except` that silently
    turns "not tested" into "absent".
    """
    if marker.derivation:
        from .derive import DERIVATIONS

        handler = DERIVATIONS.get(marker.derivation)
        if handler is None:
            raise ValueError(
                f"marker {marker.id!r} names derivation {marker.derivation!r}, "
                f"which is not registered in derive.DERIVATIONS"
            )
        return handler(genome, marker)

    rsid = marker.rsids[0]
    call = genome.get(rsid)

    if call is None:
        return Finding(marker=marker, status=FindingStatus.NOT_ON_CHIP)

    if call.is_no_call:
        return Finding(marker=marker, status=FindingStatus.NO_CALL)

    return _match(marker, call)


def _match(marker: Marker, call: Genotype) -> Finding:
    """Match a call against the catalog, flipping strand if that is what it takes.

    Consumer providers report on the plus strand of the reference, but the
    papers a marker is drawn from often report on whichever strand the original
    assay used. So a catalog entry may be keyed on `CT` while the user's file
    says `AG` — the same call, written from the other side.

    Trying the direct match first and the complement only on failure is
    load-bearing. Flipping unconditionally would corrupt every marker that was
    already correct, and the asymmetry is what makes the fallback safe: a
    genotype that matches the catalog as written is by definition already on
    the catalog's strand.

    This is safe for the A/G and C/T markers that make up most of the catalog,
    and *not* safe in general. At an A/T or C/G site the complement of one
    valid genotype is another valid genotype, so a flip cannot be detected from
    the call alone — `_is_strand_ambiguous` refuses those rather than guessing.
    """
    direct = call.canonical()
    if direct in marker.interpretations:
        return _build(marker, call, direct, flipped=False)

    if _is_strand_ambiguous(marker):
        # Both strands are internally consistent here, so a mismatch cannot be
        # resolved by flipping — we would be choosing between two answers with
        # no evidence. Report it unread instead.
        return Finding(marker=marker, status=FindingStatus.UNKNOWN_GENOTYPE, genotype=direct)

    flipped_call = call.complement()
    flipped = flipped_call.canonical()
    if flipped in marker.interpretations:
        return _build(marker, flipped_call, flipped, flipped=True)

    # Neither strand matches. Usually a genuine third allele at a multi-allelic
    # site, or a build mismatch. Surfaced rather than dropped so it can be
    # reported as an error in the catalog rather than vanishing.
    return Finding(marker=marker, status=FindingStatus.UNKNOWN_GENOTYPE, genotype=direct)


def _is_strand_ambiguous(marker: Marker) -> bool:
    """True when this marker's alleles are their own complements (A/T or C/G).

    At such a site "AT" flips to "AT" and "CG" flips to "CG", so the strand
    cannot be recovered from the genotype. Detected from the catalog's own keys
    rather than from a declared ref/alt pair, so it stays correct even for
    entries that never declare one.
    """
    alleles = {a for genotype in marker.interpretations for a in genotype if a in "ACGT"}
    return alleles in ({"A", "T"}, {"C", "G"})


def _build(marker: Marker, call: Genotype, key: str, *, flipped: bool) -> Finding:
    interpretation: Interpretation = marker.interpretations[key]
    return Finding(
        marker=marker,
        status=FindingStatus.OK,
        genotype=key,
        zygosity=_zygosity(marker, call),
        interpretation=interpretation,
        strand_flipped=flipped,
        match_method=MatchMethod.COMPLEMENT if flipped else MatchMethod.DIRECT,
    )


def _zygosity(marker: Marker, call: Genotype) -> Zygosity:
    """Count copies of the effect allele.

    `call` is expected to already be on the catalog's strand — `_match` flips it
    before calling here, so `effect_allele` compares directly.
    """
    if marker.effect_allele is None:
        return Zygosity.UNKNOWN

    alleles = call.canonical()
    if call.is_haploid:
        return Zygosity.HEMIZYGOUS if alleles == marker.effect_allele else Zygosity.HOMOZYGOUS_REF

    count = alleles.count(marker.effect_allele)
    return {
        0: Zygosity.HOMOZYGOUS_REF,
        1: Zygosity.HETEROZYGOUS,
        2: Zygosity.HOMOZYGOUS_ALT,
    }.get(count, Zygosity.UNKNOWN)
