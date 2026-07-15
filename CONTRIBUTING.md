# Contributing to DNAdistiller

The most valuable contribution is a marker: one locus, its interpretations, and the citations behind them. That is a YAML file, not code, so reviewing it needs a person who can read a paper rather than a person who can read Python.

The second most valuable is telling us a marker is wrong. Please open an issue if you find one.

## Ground rules

**Never commit a DNA file.** Not yours, not a relative's, not one from a public dataset. `tests/fixtures/README.md` explains why the rule has no exceptions. `.gitignore` blocks the obvious patterns, which stops accidents and not determination.

**Never add a network call.** The tool's central claim is that your data cannot leave the machine because there is no code that could send it. One `requests` import ends that, and it ends it for every user, not just the one who wanted the feature. A pull request that adds a network dependency will be declined regardless of what it does.

**Cite or downgrade.** A marker graded `strong` or `moderate` needs a PMID or DOI. The loader enforces this and will refuse to start otherwise. Markers graded `weak` may go uncited, since some exist specifically to record that the evidence is thin.

## Setup

```sh
git clone https://github.com/Lucaks3/DNAdistiller
cd dnadistiller
uv venv && uv pip install -e ".[dev]"

uv run pytest
uv run ruff check src tests
uv run mypy
```

## Scope: does this marker belong at all?

Answer this before anything else. It is the question most marker pull requests fail, and it fails them for good reasons.

**DNAdistiller is a longevity and healthspan tool.** A marker earns a place only if it plausibly bears on how long you stay healthy. The test:

> Would knowing this change how a reasonable person thinks about the next thirty years of their health?

That is a high bar and it is meant to be. Things that pass: cardiovascular and metabolic risk, cognitive ageing, drug safety over a lifetime, exposures that compound over decades, and the small number of variants genuinely associated with exceptional longevity.

Things that do **not** pass, however well-evidenced and however fun:

- **Traits.** Earwax type, cilantro aversion, photic sneeze, eye colour. Real genetics, correctly called, and nothing to do with healthspan.
- **Behavioural and personality genotypes.** COMT "warrior/worrier", DRD2 reward-seeking, OXTR "empathy", BDNF "plasticity". Mostly candidate-gene-era findings with poor replication records, and the ones that replicate have effects far too small to say anything about a person. Another tool in this space renders BDNF Met/Met under a red "Critical Findings" heading next to fatal myelosuppression. Do not help us become that.
- **Athletic performance.** ACTN3 is the standing example, and it is a real one: it shipped in this catalog until the scope rule above was written down, at which point it plainly failed its own test and came out. The group-level association is genuine, and it explains under 1% of sprint-time variance. A 22-author *BJSM* consensus (PMID 26582191) concluded no young athlete should be genotyped to guide training. Interesting; not healthspan.
- **Ancestry.** Different tool, different privacy model, and it needs the whole file.
- **Rare clinically actionable variants.** BRCA and similar are excluded on purpose. Consumer arrays have false-positive rates at rare positions high enough that a single call is untrustworthy, and the failure mode is somebody making an irreversible decision on a bad result from a spit tube. That needs a clinical-grade test and a genetic counsellor.

If a marker is fascinating but fails the test, the honest answer is that it belongs in a different tool. The catalog is short because the bar is high, and the bar being high is the product.

One deliberate exception exists: **MTHFR**, graded `weak`, is in the catalog on no scientific merit whatsoever. It is there as a rebuttal, because a user who arrives having been told they have an "MTHFR mutation" deserves to read what the evidence actually says. If you want to add a marker on those grounds, say so explicitly in the pull request and expect to argue for it.

## Adding a marker

Markers live in `data/markers/<category>.yaml` as a list. A minimal one:

```yaml
- id: cyp1a2_caffeine
  gene: CYP1A2
  name: Caffeine metabolism rate
  category: caffeine
  evidence: moderate
  rsids: [rs762551]
  effect_allele: C
  effect_size: "Slow metabolisers clear caffeine roughly 4x slower"
  ancestry_note: "Established mainly in European and Costa Rican cohorts."
  citations:
    - "PMID:16522833"
  interpretations:
    AA: "Fast caffeine metabolism."
    AC:
      summary: "Intermediate caffeine metabolism."
      detail: "One copy of the slow-metabolising C allele."
      population_frequency: 0.45
    CC: "Slow caffeine metabolism."
```

### Fields

| Field | Required | Notes |
|---|---|---|
| `id` | yes | Unique, snake_case. |
| `gene` | yes | Standard HGNC symbol. |
| `name` | yes | The trait in plain language, not the gene. |
| `category` | yes | One of the categories in `models.Category`. |
| `evidence` | yes | `strong`, `moderate`, or `weak`. See below. |
| `rsids` | yes | Exact dbSNP rsIDs, or 23andMe `i` ids. Verify each one. |
| `interpretations` | yes | Genotype to meaning. A string, or a mapping with `summary`, `detail`, `population_frequency`. |
| `effect_allele` | no | The allele the effect is attributed to. Drives the "one copy" phrasing. |
| `effect_size` | no | Print the real number. `OR 1.08` is more honest than "increased risk". |
| `citations` | conditional | Required for `strong` and `moderate`. |
| `ancestry_note` | no | Which populations the evidence came from. |
| `minus_strand` | no | See strand orientation below. |
| `sensitive` | no | Excluded unless `--include-sensitive`. |
| `derivation` | no | Names a handler in `derive.py`. Only for markers that are not a genotype lookup. |

### Evidence grades

`strong` means a large replicated GWAS or meta-analysis, or pharmacogenomics with actual dosing guidelines behind it. APOE and Alzheimer's qualifies. SLCO1B1 and statin myopathy qualifies.

`moderate` means replicated across independent cohorts with an effect size too small to say much about one person. Most of the catalog is here.

`weak` means small, inconsistent, or never replicated. MTHFR C677T is the reference case: enormously popular online, and the evidence for the health claims made about it does not hold up. It stays in the catalog because leaving it out reads as an oversight, and because a user who has been told online that they have a "MTHFR mutation" deserves to see what the evidence actually says.

If you are unsure between two grades, take the lower one. A tool that undersells a real finding costs someone a conversation with their GP. A tool that oversells a weak one costs someone a supplement habit and a false belief about their body, and it costs the project the only thing it has, which is being trustworthy about uncertainty.

### Strand orientation

This is the most common way to get a marker wrong, and it fails silently.

Consumer providers report on the plus strand of the reference genome. Papers report on whichever strand the original assay used. So a paper may describe a variant as `C/T` while a user's file says `G/A` for the same call.

Write your `interpretations` keys on the plus strand, matching what a 23andMe file would contain. If your source reports on the minus strand, flip the alleles yourself (A to T, C to G) and set `minus_strand: true` to record that you did. Do not transcribe the paper's letters directly and hope.

The interpreter tries a direct match first and only then the complement, so an A/G or C/T marker written on the wrong strand will usually still resolve. It cannot save you at an A/T or C/G site, where the complement of a valid genotype is another valid genotype. At those sites the tool refuses to guess and reports the marker unread, so getting the strand right in the YAML is the only fix.

Check the plus-strand alleles at [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) before you submit. It takes a minute and it is the difference between a result and a wrong result.

### The inversion, which is the one that will get you

Getting the strand wrong is recoverable. Getting the alleles the wrong way round is not, and it is the failure mode that has actually hurt people.

A published consumer-DNA tool shipped TPMT rs1142345 with `TT` labelled *"poor metaboliser, standard thiopurine doses can cause fatal myelosuppression"*. `T` is the reference allele there, carried by about 96% of people. The entry was backwards, so nearly every user was told they risked a fatal drug reaction. It reached their sample report, and their own test suite asserted the wrong result, which is what let it survive.

The counterintuitive part is that a strand fallback makes this **worse**, not better. Complement matching means your orientation choice no longer affects whether a genotype *matches*, while it still entirely decides *what the match means*. A minus-strand file complements neatly onto an inverted entry and gets the same wrong answer. Strand tolerance hides orientation errors instead of surfacing them.

You cannot catch this by reading. An inverted entry looks completely reasonable unless you already know which allele is the reference. So the validator catches it mechanically, by comparing your `effect_allele` against your declared `population_frequency`: an effect allele whose homozygote is carried by most of the population is almost always inverted.

That check only works if you fill in `population_frequency`. Please do. It is the field that makes the worst bug in this domain unrepresentable rather than merely unlikely.

### Verify by retrieval, never from memory

Every citation must be checked by actually opening it. Not recalled, not inferred from the title.

This is not pedantry. During this catalog's construction, three of three PMIDs written from memory turned out to be wrong, and two of them resolved to real, plausible-looking papers on adjacent topics. A wrong PMID that resolves to a real paper is worse than no citation: it survives review, because the reviewer clicks it and sees a paper.

### Writing interpretations

Write `summary` so it stands alone. At minimal tier it is all the reader gets, stripped of the gene name and the rsID, so "Reduced enzyme activity" works and "This variant reduces it" does not.

Describe, do not prescribe. "Slow caffeine metabolism" is a finding. "You should avoid coffee after 2pm" is advice, and this tool does not give advice. The LLM the profile gets pasted into can have that conversation with the user, with their context, which is the whole design.

Include the common genotype too. A person whose result is "typical activity" has learned something, and a catalog that only describes the risk allele produces profiles where every entry sounds like bad news.

## Adding a provider

Add a `Parser` subclass in `src/dnadistiller/parsers/providers.py` and list it in `PARSERS`. Keep `matches` narrow: these formats look alike, and a loose matcher claims a competitor's file and misreads every column in it. Add a synthetic fixture.

`PARSERS` is ordered. More specific parsers come first.

## Pull requests

Run `pytest`, `ruff check`, and `mypy` before opening one. CI runs all three.

For a marker, say in the description where the effect allele and strand came from. A link to the dbSNP page is ideal.
