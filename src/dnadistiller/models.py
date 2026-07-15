"""Core domain types for dnadistiller.

The vocabulary here maps onto the pipeline the tool implements:

    raw export file  ->  Genome    (every call the chip made, ~600k of them)
    marker catalog   ->  Marker    (a curated locus we know how to interpret)
    Genome x Marker  ->  Finding   (what this person's calls mean at that locus)
    [Finding]        ->  Profile   (the redacted, shareable subset)

`Genome` is the only object that ever holds the full raw dataset, and it is
never serialised to disk. Everything downstream of `Finding` is deliberately
lossy: see `Tier` for the reason.
"""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field
from enum import Enum

# ---------------------------------------------------------------------------
# Alleles and genotypes
# ---------------------------------------------------------------------------

#: Watson-Crick complements, used to flip a call to the other strand.
#: I and D are 23andMe's insertion/deletion pseudo-alleles and have no
#: complement — they are strand-agnostic, so they map to themselves.
COMPLEMENT: dict[str, str] = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "I": "I",
    "D": "D",
    "-": "-",
    "0": "0",
}

#: Tokens a provider may use to mean "the chip did not produce a call here".
NO_CALL_TOKENS: frozenset[str] = frozenset({"--", "-", "00", "0", "", "NN", "N", "NC"})


class Zygosity(Enum):
    """How many copies of the non-reference allele are present.

    Reported at STANDARD tier as a stand-in for the raw genotype: "one copy"
    carries most of the interpretive value of "AG" while being a coarser,
    less identifying statement.
    """

    HOMOZYGOUS_REF = "no copies"
    HETEROZYGOUS = "one copy"
    HOMOZYGOUS_ALT = "two copies"
    HEMIZYGOUS = "one copy (single chromosome)"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class Genotype:
    """A single call from the array.

    `alleles` is stored exactly as the provider wrote it, uppercased. It is not
    normalised to a strand at parse time on purpose: the correct strand depends
    on which marker you are asking about, so normalisation happens at
    interpretation time where the catalog can say what it expects.
    """

    rsid: str
    chromosome: str
    position: int
    alleles: str

    @property
    def is_no_call(self) -> bool:
        return self.alleles.upper() in NO_CALL_TOKENS or set(self.alleles) <= {"-", "0"}

    @property
    def is_haploid(self) -> bool:
        """Single-allele call: male X/Y outside the pseudoautosomal regions, and MT."""
        return len(self.alleles) == 1 and not self.is_no_call

    @property
    def is_ambiguous(self) -> bool:
        """Neither a clean genotype nor a recognised no-call.

        Counted at parse time as a health signal rather than acted on. One of
        these is a curiosity; thousands mean a delimiter or encoding problem is
        mangling the last column, and that is worth telling the user before they
        read any results derived from it.
        """
        if self.is_no_call:
            return False
        alleles = self.alleles.upper()
        if len(alleles) not in (1, 2):
            return True
        return any(allele not in "ACGTID" for allele in alleles)

    def complement(self) -> Genotype:
        """Return this call as it would read on the opposite strand.

        Alleles are re-sorted afterwards because genotype strings are unordered
        pairs — the array does not tell us which allele sits on which
        chromosome, so "AG" and "GA" denote the same call and we canonicalise.
        """
        if self.is_no_call:
            return self
        flipped = "".join(COMPLEMENT.get(a, a) for a in self.alleles.upper())
        return Genotype(
            rsid=self.rsid,
            chromosome=self.chromosome,
            position=self.position,
            alleles="".join(sorted(flipped)),
        )

    def canonical(self) -> str:
        """Alleles sorted alphabetically, so AG and GA compare equal."""
        return "".join(sorted(self.alleles.upper()))


# ---------------------------------------------------------------------------
# Genome
# ---------------------------------------------------------------------------


class Build(Enum):
    """Reference genome assembly.

    Positions are meaningless without one of these. We key lookups on rsID
    rather than position specifically so that a build mismatch degrades into a
    missing marker rather than a wrong answer at the wrong locus. That choice is
    what lets this tool treat build as provenance to report rather than a
    coordinate system to remap, and it is why an old export still works.

    NCBI36 is here because it is real: 23andMe exports downloaded around 2011
    declare `build 36`, and the same v3 chip exported build 36 in 2011 and build
    37 in 2014. Build is a property of when the file was downloaded, not of the
    chip. Positions moved by megabases between the two (rs429358 sits at
    19:50103781 in build 36 and 19:45411941 in build 37), so anything reading
    positions must not assume 37.
    """

    NCBI36 = 36
    GRCH37 = 37
    GRCH38 = 38
    UNKNOWN = 0


class Severity(Enum):
    INFO = "info"
    WARNING = "warning"
    ERROR = "error"

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True, slots=True)
class ParseIssue:
    """Something worth telling the user about their file.

    Structured rather than printed at the point of discovery, because the same
    facts need to reach a terminal, a profile, and a future GUI, and because a
    `code` is greppable in a bug report where prose is not.

    Not exceptions: none of these stop a parse. They are the difference between
    a file that read cleanly and one that read but should be distrusted, which
    is a distinction the user cannot make for themselves.
    """

    severity: Severity
    code: str
    message: str


@dataclass
class Genome:
    """Every call in a person's raw export.

    This is the sensitive object. It stays in memory, is never written to disk
    by this tool, and nothing that leaves the process should be derived from it
    except through `Profile`.
    """

    genotypes: dict[str, Genotype] = field(default_factory=dict)
    build: Build = Build.UNKNOWN
    source: str = "unknown"

    #: Rows the parser could not read at all. A handful is normal; thousands
    #: means the delimiter or column layout was misread and every result
    #: downstream is suspect.
    malformed_rows: int = 0

    #: rsIDs appearing more than once. Real exports contain these. The last
    #: occurrence wins, which is arbitrary, so the count is surfaced rather
    #: than the choice being made silently.
    duplicate_rsids: int = 0

    #: Calls that are neither a clean genotype nor a recognised no-call token.
    #: Counted at parse time because a spike here is the signature of a
    #: delimiter or encoding problem rather than of biology.
    ambiguous_calls: int = 0

    #: Diagnostics for the user. See `ParseIssue`.
    issues: list[ParseIssue] = field(default_factory=list)

    #: What this provider's array structurally cannot see, regardless of result.
    #: Carried on the genome rather than looked up at render time so that the
    #: caveat travels with the data it qualifies.
    blind_spots: list[str] = field(default_factory=list)

    #: Provider-reported sex chromosomes present, used to sanity-check haploid calls.
    _chromosomes_seen: set[str] = field(default_factory=set, repr=False)

    def __len__(self) -> int:
        return len(self.genotypes)

    def __contains__(self, rsid: str) -> bool:
        return rsid.lower() in self.genotypes

    def __iter__(self) -> Iterator[Genotype]:
        return iter(self.genotypes.values())

    def get(self, rsid: str) -> Genotype | None:
        """Look up a call by rsID, case-insensitively.

        Returns None when the variant is absent from this chip version, which
        is a routine and expected outcome rather than an error: chip content
        varies by provider and by year, and no export contains every marker in
        the catalog.
        """
        return self.genotypes.get(rsid.lower())

    def add(self, genotype: Genotype) -> None:
        """Record a call, counting the two things worth knowing about it.

        Duplicate rsIDs are real in real exports, and last-one-wins is an
        arbitrary choice. It is made here once, and counted, so the arbitrariness
        is reported rather than hidden.
        """
        key = genotype.rsid.lower()
        if key in self.genotypes:
            self.duplicate_rsids += 1
        if genotype.is_ambiguous:
            self.ambiguous_calls += 1
        self.genotypes[key] = genotype
        self._chromosomes_seen.add(genotype.chromosome)

    @property
    def has_y_chromosome(self) -> bool:
        return "Y" in self._chromosomes_seen

    @property
    def chromosomes(self) -> set[str]:
        return set(self._chromosomes_seen)


# ---------------------------------------------------------------------------
# Marker catalog
# ---------------------------------------------------------------------------


class Evidence(Enum):
    """How much weight a marker's claim can bear.

    Consumer genomics is full of confidently-stated findings resting on one
    small study that never replicated. Grading every marker forces that
    judgement to be explicit and reviewable in a pull request, and lets users
    filter the catalog down to what actually holds up.
    """

    STRONG = "strong"
    """Large replicated GWAS or meta-analysis, or clinically actionable
    pharmacogenomics with dosing guidelines behind it."""

    MODERATE = "moderate"
    """Replicated across independent cohorts, but with a modest effect size
    that says little about any individual."""

    WEAK = "weak"
    """Small, inconsistent, or unreplicated. Included only where a marker is
    popular enough that its absence would be read as an oversight — shown with
    the evidence against it stated plainly."""

    def __str__(self) -> str:
        return self.value


class Category(Enum):
    """Topic grouping, and the unit users opt in and out of."""

    LIPIDS = "lipids"
    METABOLIC = "metabolic"
    CLOTTING = "clotting"
    RESPIRATORY = "respiratory"
    METHYLATION = "methylation"
    CAFFEINE = "caffeine"
    INFLAMMATION = "inflammation"
    NEURO = "neuro"
    LONGEVITY = "longevity"
    FITNESS = "fitness"
    SLEEP = "sleep"
    PHARMACOGENOMICS = "pharmacogenomics"

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True, slots=True)
class Interpretation:
    """What one genotype at one marker means.

    `summary` is written to survive being stripped of all identifiers at
    MINIMAL tier, so it has to stand alone as a sentence: "Slow caffeine
    metaboliser" rather than "This variant slows it down".
    """

    summary: str
    detail: str = ""
    #: Frequency of this genotype in the general population, 0..1, where known.
    #: Present so a result can be framed as common or unusual rather than as
    #: an unqualified verdict.
    population_frequency: float | None = None


class TestVerdict(Enum):
    """How a measurement relates to the genotype that suggested it.

    The three-way split is the honest core of this tool. DNA is bad at "how am I
    doing" and good at "what is worth measuring", and being explicit about which
    of these three a marker is stops the catalog from quietly overselling itself.
    """

    SUPERSEDES = "supersedes"
    """The test beats the genotype outright. Ferritin over HFE: iron loading is
    what harms you, and most C282Y homozygotes never load iron. Most of this
    catalog is here, which is a fact about genetics rather than a shortcoming."""

    COMPLEMENTS = "complements"
    """Worth having, but it answers a different question than the variant does."""

    NONE_EXISTS = "none"
    """No test warns you in time. Nothing in a panel gives you your Factor V
    Leiden status before the clot, and CK only rises after the muscle damage.
    This is the small set where the genotype is the whole point."""

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True, slots=True)
class BloodTest:
    """The measurement a genotype is only a prior for.

    A genotype gives you a lifetime prior. A blood test gives you the actual
    number, today, for about the price of lunch. A tool that reports the prior
    and never names the measurement has stopped halfway.
    """

    verdict: TestVerdict
    #: Empty when `verdict` is NONE_EXISTS: there is nothing to name.
    name: str = ""
    note: str = ""

    def __post_init__(self) -> None:
        if self.verdict is not TestVerdict.NONE_EXISTS and not self.name:
            raise ValueError(f"a {self.verdict} blood_test must name the test")
        if self.verdict is TestVerdict.NONE_EXISTS and self.name:
            raise ValueError(
                "a blood_test with verdict 'none' must not name a test; "
                "put the explanation in `note`"
            )


@dataclass(frozen=True, slots=True)
class Marker:
    """One curated locus, loaded from `data/markers/*.yaml`."""

    id: str
    gene: str
    name: str
    category: Category
    evidence: Evidence
    rsids: tuple[str, ...]
    interpretations: dict[str, Interpretation]

    summary: str = ""
    effect_size: str = ""
    citations: tuple[str, ...] = ()

    #: The allele the literature attributes the effect to, on the same strand as
    #: `interpretations`. Counting its copies is what turns a genotype into the
    #: "one copy / two copies" phrasing that STANDARD tier reports in place of
    #: the raw call. Left unset for markers where "effect allele" does not apply
    #: — APOE's haplotypes being the obvious one.
    effect_allele: str | None = None

    #: Strand the catalog's genotype keys are written on, relative to the
    #: reference plus strand. When a source paper reports on the minus strand we
    #: record that here and flip the user's call before matching, rather than
    #: silently transcribing flipped keys into the YAML where nobody can audit
    #: the decision later.
    minus_strand: bool = False

    #: Findings that carry weight beyond curiosity — APOE and Alzheimer's being
    #: the canonical case. Excluded unless explicitly requested.
    sensitive: bool = False

    #: Free-text note on which ancestries the evidence was established in.
    #: Most GWAS is European-ancestry-biased and transfers poorly; saying so
    #: per-marker is more honest than one blanket disclaimer nobody reads.
    ancestry_note: str = ""

    #: Set when a marker needs logic rather than a genotype lookup, e.g. APOE's
    #: e2/e3/e4 haplotype derived from two SNPs. Names a handler in derive.py.
    derivation: str | None = None

    #: The measurement worth having instead of, or alongside, this genotype.
    #: See `BloodTest`. Optional, but a marker without one should have a reason:
    #: usually that no test exists until the harm has already happened.
    blood_test: BloodTest | None = None

    def __post_init__(self) -> None:
        if not self.rsids:
            raise ValueError(f"marker {self.id!r} declares no rsIDs")
        if not self.derivation and len(self.rsids) != 1:
            raise ValueError(
                f"marker {self.id!r} lists {len(self.rsids)} rsIDs but no derivation; "
                "multi-SNP markers need a handler to combine them"
            )


# ---------------------------------------------------------------------------
# Findings and profile
# ---------------------------------------------------------------------------


class FindingStatus(Enum):
    OK = "ok"
    NOT_ON_CHIP = "not_on_chip"
    """The provider's array never tested this position."""
    NO_CALL = "no_call"
    """Tested, but the assay failed to produce a call."""
    UNKNOWN_GENOTYPE = "unknown_genotype"
    """Called, but the result is not one the catalog knows how to read —
    usually a strand or build problem, so it is surfaced rather than dropped."""


class MatchMethod(Enum):
    """How a call was matched to a catalog entry.

    Recorded rather than discarded because the methods are not equally certain,
    and the difference is invisible in the result. A direct match is an
    observation. A complement match is an *inference* that the file reported on
    the other strand, and while that inference is sound at a non-palindromic
    site, it is still a step of reasoning the reader deserves to see. Silently
    complementing and presenting both as the same kind of fact is how a tool
    ends up confidently wrong without anyone able to tell.
    """

    DIRECT = "direct"
    COMPLEMENT = "complement"
    DERIVED = "derived"
    """Computed from several SNPs by a handler in derive.py, e.g. APOE."""

    def __str__(self) -> str:
        return self.value


@dataclass(frozen=True, slots=True)
class Finding:
    marker: Marker
    status: FindingStatus
    genotype: str | None = None
    zygosity: Zygosity = Zygosity.UNKNOWN
    interpretation: Interpretation | None = None
    #: Set when the call had to be flipped to match the catalog's strand.
    #: Surfaced in verbose output so the flip is inspectable.
    strand_flipped: bool = False
    #: How this finding was arrived at. See `MatchMethod`.
    match_method: MatchMethod = MatchMethod.DIRECT

    @property
    def is_reportable(self) -> bool:
        return self.status is FindingStatus.OK and self.interpretation is not None


class Tier(Enum):
    """How much detail a profile discloses.

    The tiers exist because "share your genome with an LLM" and "share nothing"
    are not the only options, and the middle ground is where this tool lives.
    Each step down strips identifying detail while keeping most of what makes
    the profile worth discussing:

      FULL     rs1801133  CT  MTHFR   Reduced enzyme activity
      STANDARD MTHFR, one copy       Reduced enzyme activity
      MINIMAL  Reduced MTHFR enzyme activity

    None of these are anonymous — see docs/privacy.md. Genotypes at even a few
    dozen independent common SNPs are enough to single someone out, so the
    honest claim is minimisation, not anonymisation. What the tiers buy you is
    control over how much you hand over, and a floor on the damage if the
    recipient turns out to be less trustworthy than you assumed.
    """

    MINIMAL = "minimal"
    STANDARD = "standard"
    FULL = "full"

    def __str__(self) -> str:
        return self.value


@dataclass
class Profile:
    """The redacted, shareable output. The only object intended to leave the machine."""

    findings: list[Finding]
    tier: Tier
    source: str
    build: Build
    #: Markers requested but absent from this chip, kept so the profile can say
    #: "not tested" instead of leaving a silent gap that reads as "no risk".
    missing: list[Marker] = field(default_factory=list)

    #: What the source array structurally cannot see. Travels into the shared
    #: profile at every tier, because it is a fact about the chip rather than
    #: about the person: it discloses nothing and it is the context most likely
    #: to stop a reader treating silence as an all-clear.
    blind_spots: list[str] = field(default_factory=list)

    #: Parse diagnostics worth passing on, notably a truncated source file. A
    #: reader who does not know the input was 4,000 variants instead of 600,000
    #: will read the gaps as biology.
    issues: list[ParseIssue] = field(default_factory=list)

    @property
    def reportable(self) -> list[Finding]:
        return [f for f in self.findings if f.is_reportable]
