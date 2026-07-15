"""The preamble and follow-up prompts that ship alongside a profile.

A profile pasted into a chatbot with no framing tends to get read badly. The
model sees a list of risk alleles and obliges with a confident narrative, and an
odds ratio of 1.12 comes back as "you are at elevated risk for X" — which is not
what the data says, and is the specific failure this module exists to prevent.

So the profile leads with a preamble aimed at the model rather than the user. It
is not a disclaimer for our benefit; it is context that measurably changes the
answer, and it is why `dnadistiller profile` emits it by default instead of leaving
users to think of it themselves.
"""

from __future__ import annotations

from .models import Category, Profile, Tier

PREAMBLE = """\
This is a partial genetic profile, generated locally by an open-source tool
(DNAdistiller) from a consumer DNA test. Please read the context below before
interpreting it. It changes what these results can support.

What this data is:
- A hand-picked subset of common variants related to longevity and healthspan.
  It is not a genome, not an exome, and not a clinical test.
- Output from a consumer genotyping array. Arrays test preselected positions.
  They do not sequence. Anything absent from this profile was probably never
  tested, and "not tested" does not mean "normal".
- Not clinically validated. At rare positions these arrays have a false-positive
  rate high enough that a single finding is not trustworthy on its own.

What this data is not:
- A basis for predictions about this person. Most of these variants have odds
  ratios between 1.05 and 1.3. At that size a variant says almost nothing about
  one individual, even where the association is real and replicated.
- Causal. Association is not causation, and neither is diagnosis.
- Necessarily applicable. Most of the underlying research used
  European-ancestry cohorts and often transfers poorly to other groups.
- The main driver. Genetics contributes modestly to most longevity outcomes
  next to smoking, sleep, exercise, diet, alcohol, and access to healthcare.

How to be useful here:
- Be calibrated rather than reassuring or alarming. If a finding is weak, say
  it is weak. "This probably does not matter much for you" is a useful answer,
  not a failure to engage.
- Separate what is actionable from what is merely interesting. Very little of
  this is actionable.
- Where something is worth following up, name who with: a GP, a genetic
  counsellor, a lipid clinic. Do not recommend supplements or dosing.
- On medication, defer to a clinician. Pharmacogenomic markers can inform a
  prescribing conversation. They cannot replace one.
"""

TIER_NOTES: dict[Tier, str] = {
    Tier.MINIMAL: (
        "Disclosure tier: MINIMAL. Interpretations only, with no rsIDs, genotypes, or "
        "gene names. You cannot look these variants up. That is intentional. Work from "
        "the summaries as written."
    ),
    Tier.STANDARD: (
        "Disclosure tier: STANDARD. Gene names and copy counts, no raw genotypes or rsIDs."
    ),
    Tier.FULL: (
        "Disclosure tier: FULL. rsIDs and genotypes are included, so findings can be "
        "checked against the literature."
    ),
}

#: Follow-ups worth asking regardless of what the profile contains.
_GENERAL_PROMPTS: list[str] = [
    "Which of these findings, if any, would actually change something I do? "
    "Be strict. I would rather hear 'none of them' than a list of maybes.",
    "Rank these by how much evidence is really behind them. Where am I at risk of "
    "over-reading a weak result?",
    "For each finding, roughly how common is this genotype? I want to know whether "
    "I am unusual or ordinary here.",
    "What lifestyle factors have a bigger effect on my healthspan than anything in "
    "this profile? I want to keep this in proportion.",
    "Are any of these worth confirming with a clinical-grade test rather than a consumer array?",
]

#: Follow-ups that only make sense when a category is actually present.
_CATEGORY_PROMPTS: dict[Category, list[str]] = {
    Category.LIPIDS: [
        "How do these lipid-related variants relate to what a standard blood panel "
        "would show? Would a lipid panel or Lp(a) test tell me more than this does?",
        "Is there anything here that would justify testing Lp(a) once, given it is "
        "largely genetically determined and rarely measured?",
    ],
    Category.METABOLIC: [
        "Do these metabolic variants suggest anything beyond the standard advice on "
        "diet and exercise, or are they just describing baseline risk?",
        "How would these compare against an actual HbA1c or fasting glucose result?",
    ],
    Category.CAFFEINE: [
        "Given these caffeine metabolism variants, what would a sensible experiment "
        "look like: timing, dose, and what I would measure?",
        "Does my caffeine metabolism interact with sleep quality in a way worth testing?",
    ],
    Category.METHYLATION: [
        "I have seen a lot of strong claims about MTHFR online. What does the evidence "
        "actually support, and what is overstated?",
        "Would a homocysteine or B12 blood test tell me more than the genotype does?",
    ],
    Category.PHARMACOGENOMICS: [
        "Which of these pharmacogenomic findings would be worth mentioning to a doctor "
        "before a prescription, and at what point in the conversation?",
        "Are there CPIC or equivalent guidelines covering these variants? What do they say?",
    ],
    Category.NEURO: [
        "Help me put this in proportion against age, family history, and lifestyle. "
        "What is the absolute risk here, not the relative risk?",
        "What is actually supported for cognitive healthspan, independent of genotype?",
    ],
    Category.FITNESS: [
        "Do these exercise-related variants justify changing how I train, or is the "
        "effect too small to act on?",
    ],
    Category.SLEEP: [
        "How do these circadian variants square with my actual sleep patterns? What "
        "would be worth tracking?",
    ],
    Category.INFLAMMATION: [
        "How do these inflammation variants compare with simply measuring hs-CRP?",
    ],
    Category.LONGEVITY: [
        "These are described as longevity-associated. How strong is that evidence "
        "really, and what does it mean for an individual rather than a population?",
    ],
}


def suggest_prompts(profile: Profile, *, limit: int = 8) -> list[str]:
    """Build follow-up prompts tailored to what a profile actually contains.

    Category prompts come first: they are the ones a user would not have thought
    to ask, and they are the reason this is generated per-profile rather than
    printed from a fixed list. General prompts backfill to `limit`.

    At MINIMAL tier the category prompts are dropped entirely, because several of
    them name the gene they are about ("What does the evidence on MTHFR actually
    support?"). A profile that carefully withholds the gene name and then asks a
    question containing it has disclosed the gene name. The questions are part of
    the payload, so they redact with everything else.
    """
    prompts: list[str] = []

    if profile.tier is not Tier.MINIMAL:
        present = {f.marker.category for f in profile.reportable}
        for category in sorted(present, key=str):
            prompts.extend(_CATEGORY_PROMPTS.get(category, []))

    for prompt in _GENERAL_PROMPTS:
        if len(prompts) >= limit:
            break
        prompts.append(prompt)

    return prompts[:limit]


def build_preamble(profile: Profile) -> str:
    """Preamble plus the tier note describing how redacted this profile is."""
    return f"{PREAMBLE}\n{TIER_NOTES[profile.tier]}\n"
