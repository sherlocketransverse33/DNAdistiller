"""Loading and validating the marker catalog.

The catalog lives in `data/markers/*.yaml`, not in Python. That split is the
main thing this package exists to enforce, and it is deliberate:

- A marker is a claim about human biology backed by citations. Reviewing one
  means checking an rsID against dbSNP and a claim against a paper — a job for
  a domain reader, who should not have to read Python to do it.
- Contributors adding a marker should be writing data, not code. A pull request
  that is ten lines of YAML with two PMIDs is one a geneticist can review on
  their phone.
- The validation here is strict for the same reason. A typo in an rsID does not
  fail loudly at runtime; it silently reports "not tested" forever. So the
  errors below are pedantic on purpose, and they name the file and marker,
  because they are read by contributors rather than by us.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from ..models import BloodTest, Category, Evidence, Interpretation, Marker, TestVerdict

#: Bundled with the wheel via pyproject's force-include; falls back to the repo
#: layout when running from a source checkout.
_PACKAGED = Path(__file__).parent.parent / "_data" / "markers"
_REPO = Path(__file__).parent.parent.parent.parent / "data" / "markers"


class CatalogError(Exception):
    """Raised when a catalog file is malformed. Message names the offending marker."""


def catalog_dir() -> Path:
    if _PACKAGED.is_dir():
        return _PACKAGED
    if _REPO.is_dir():
        return _REPO
    raise CatalogError(
        f"no marker catalog found at {_PACKAGED} or {_REPO}. "
        "If running from a source checkout, ensure data/markers/ exists."
    )


def load_catalog(path: Path | None = None) -> list[Marker]:
    """Load every marker from the catalog directory.

    Markers are returned sorted by id so that output ordering is stable across
    runs and filesystems — a profile that reshuffles between runs is annoying to
    diff and makes the tool look nondeterministic.
    """
    directory = path or catalog_dir()
    markers: list[Marker] = []
    seen: dict[str, Path] = {}

    for yaml_file in sorted(directory.glob("*.yaml")):
        for raw in _read_file(yaml_file):
            marker = _parse_marker(raw, yaml_file)
            if marker.id in seen:
                raise CatalogError(
                    f"duplicate marker id {marker.id!r} in {yaml_file.name}; "
                    f"already defined in {seen[marker.id].name}"
                )
            seen[marker.id] = yaml_file
            markers.append(marker)

    if not markers:
        raise CatalogError(f"catalog directory {directory} contains no markers")

    return sorted(markers, key=lambda m: m.id)


def _read_file(path: Path) -> list[dict[str, Any]]:
    try:
        # safe_load, not load: the catalog is data from pull requests, and
        # yaml.load on untrusted input constructs arbitrary Python objects.
        content = yaml.safe_load(path.read_text(encoding="utf-8"))
    except yaml.YAMLError as exc:
        raise CatalogError(f"{path.name} is not valid YAML: {exc}") from exc

    if content is None:
        return []
    if not isinstance(content, list):
        raise CatalogError(
            f"{path.name} must contain a list of markers at the top level, "
            f"got {type(content).__name__}"
        )
    return content


def _parse_marker(raw: dict[str, Any], source: Path) -> Marker:
    if not isinstance(raw, dict):
        raise CatalogError(f"{source.name}: expected a marker mapping, got {type(raw).__name__}")

    marker_id = raw.get("id")
    if not marker_id:
        raise CatalogError(f"{source.name}: a marker is missing its required 'id' field")

    where = f"{source.name}:{marker_id}"

    for required in ("gene", "name", "category", "evidence", "rsids", "interpretations"):
        if required not in raw:
            raise CatalogError(f"{where}: missing required field {required!r}")

    try:
        category = Category(raw["category"])
    except ValueError:
        valid = ", ".join(sorted(c.value for c in Category))
        raise CatalogError(
            f"{where}: unknown category {raw['category']!r}. Valid categories: {valid}"
        ) from None

    try:
        evidence = Evidence(raw["evidence"])
    except ValueError:
        valid = ", ".join(sorted(e.value for e in Evidence))
        raise CatalogError(
            f"{where}: unknown evidence grade {raw['evidence']!r}. Valid grades: {valid}"
        ) from None

    rsids = raw["rsids"]
    if isinstance(rsids, str):
        rsids = [rsids]
    if not isinstance(rsids, list) or not rsids:
        raise CatalogError(f"{where}: 'rsids' must be a non-empty list")
    for rsid in rsids:
        _validate_rsid(rsid, where)

    citations = raw.get("citations") or []
    if isinstance(citations, str):
        citations = [citations]
    if evidence in (Evidence.STRONG, Evidence.MODERATE) and not citations:
        # Weak markers are allowed to be uncited because some are included
        # precisely to document that the evidence is thin. A strong claim with
        # no citation is how folklore gets laundered into a catalog.
        raise CatalogError(
            f"{where}: evidence is {evidence.value!r} but no citations are given. "
            "Markers graded strong or moderate must cite a PMID or DOI."
        )

    # Derived markers key their interpretations on names their handler produces
    # ("c282y_homozygous", "e3/e4") rather than on allele pairs, so the allele
    # validation below does not apply to them.
    interpretations = _parse_interpretations(
        raw["interpretations"], where, derived=bool(raw.get("derivation"))
    )

    effect_allele = raw.get("effect_allele")
    if effect_allele is not None:
        effect_allele = str(effect_allele).upper()
        if effect_allele not in {"A", "C", "G", "T", "I", "D"}:
            raise CatalogError(
                f"{where}: effect_allele {effect_allele!r} is not a single valid allele"
            )

    blood_test = _parse_blood_test(raw.get("blood_test"), where)

    try:
        return Marker(
            id=str(marker_id),
            gene=str(raw["gene"]),
            name=str(raw["name"]),
            category=category,
            evidence=evidence,
            rsids=tuple(str(r) for r in rsids),
            interpretations=interpretations,
            summary=str(raw.get("summary", "")),
            effect_size=str(raw.get("effect_size", "")),
            citations=tuple(str(c) for c in citations),
            effect_allele=effect_allele,
            minus_strand=bool(raw.get("minus_strand", False)),
            sensitive=bool(raw.get("sensitive", False)),
            ancestry_note=str(raw.get("ancestry_note", "")),
            derivation=raw.get("derivation"),
            blood_test=blood_test,
        )
    except ValueError as exc:
        # Surfaces Marker.__post_init__'s structural checks with file context.
        raise CatalogError(f"{where}: {exc}") from exc


def _parse_blood_test(raw: Any, where: str) -> BloodTest | None:
    """Parse the `blood_test:` block.

    A bare string is allowed for the common case where the test name says it all
    and it does not supersede the genotype.
    """
    if raw is None:
        return None
    if isinstance(raw, str):
        return BloodTest(verdict=TestVerdict.COMPLEMENTS, name=raw)
    if not isinstance(raw, dict):
        raise CatalogError(f"{where}: 'blood_test' must be a string or a mapping")

    try:
        verdict = TestVerdict(raw.get("verdict", "complements"))
    except ValueError:
        valid = ", ".join(v.value for v in TestVerdict)
        raise CatalogError(
            f"{where}: unknown blood_test verdict {raw.get('verdict')!r}. Valid: {valid}"
        ) from None

    try:
        return BloodTest(
            verdict=verdict,
            name=str(raw.get("name", "")),
            note=str(raw.get("note", "")),
        )
    except ValueError as exc:
        raise CatalogError(f"{where}: {exc}") from exc


def _validate_rsid(rsid: Any, where: str) -> None:
    """Check an identifier looks like an rsID or a 23andMe internal ID.

    23andMe assigns `i`-prefixed ids to variants without an rsID on the chip, so
    both forms are legitimate. A typo caught here is worth a great deal: an
    invalid rsID is invisible at runtime, silently reporting 'not tested' for
    every user forever.
    """
    if not isinstance(rsid, str):
        raise CatalogError(f"{where}: rsid {rsid!r} must be a string")
    lowered = rsid.lower()
    if lowered.startswith("rs"):
        if not lowered[2:].isdigit():
            raise CatalogError(
                f"{where}: {rsid!r} is not a valid rsID (expected rs followed by digits)"
            )
    elif lowered.startswith("i"):
        if not lowered[1:].isdigit():
            raise CatalogError(f"{where}: {rsid!r} is not a valid 23andMe internal id")
    else:
        raise CatalogError(
            f"{where}: {rsid!r} is not a recognised identifier "
            "(expected an rsID like rs1801133, or a 23andMe id like i3003137)"
        )


def _parse_interpretations(
    raw: Any, where: str, *, derived: bool = False
) -> dict[str, Interpretation]:
    if not isinstance(raw, dict) or not raw:
        raise CatalogError(f"{where}: 'interpretations' must be a non-empty mapping")

    parsed: dict[str, Interpretation] = {}
    for genotype, value in raw.items():
        key = str(genotype).strip() if derived else _canonical_genotype(str(genotype), where)

        if isinstance(value, str):
            parsed[key] = Interpretation(summary=value)
            continue

        if not isinstance(value, dict):
            raise CatalogError(
                f"{where}: interpretation for {genotype!r} must be a string or a mapping"
            )
        if "summary" not in value:
            raise CatalogError(f"{where}: interpretation for {genotype!r} needs a 'summary'")

        frequency = value.get("population_frequency")
        if frequency is not None:
            frequency = float(frequency)
            if not 0.0 <= frequency <= 1.0:
                raise CatalogError(
                    f"{where}: population_frequency for {genotype!r} is {frequency}, "
                    "expected a proportion between 0 and 1"
                )

        parsed[key] = Interpretation(
            summary=str(value["summary"]),
            detail=str(value.get("detail", "")),
            population_frequency=frequency,
        )

    return parsed


def _canonical_genotype(genotype: str, where: str) -> str:
    """Normalise a catalog genotype key to sorted uppercase.

    Genotype strings from an array are unordered pairs — the chip does not
    resolve phase — so `AG` and `GA` name the same call. Sorting both the
    catalog keys here and the user's calls at lookup time means the two always
    meet in the same form, and a contributor writing `GA` gets the same
    behaviour as one writing `AG` instead of a silent miss.
    """
    cleaned = genotype.strip().upper()

    for allele in cleaned:
        if allele not in "ACGTID-0":
            raise CatalogError(
                f"{where}: genotype key {genotype!r} contains {allele!r}, "
                "which is not a valid allele"
            )
    return "".join(sorted(cleaned))
