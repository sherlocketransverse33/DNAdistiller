# DNAdistiller

**Distil your genome down to what matters, and share only that.**

Turn a 23andMe or AncestryDNA file into a short, longevity-focused profile you can actually discuss with an LLM, without handing over your entire genome. Your raw file has 640,000 variants. A distillate has twelve.

Runs locally. No account, no upload, no network code.

```console
$ dnadistiller profile ~/genome.txt --tier standard

# Longevity genetic profile
...
## Findings

### Lipids
- **LPA** (Lipoprotein(a) level): No copies of the Lp(a)-raising variant at this locus. (no copies, evidence: strong)

### Longevity
- **FOXO3** (Longevity-associated variant): One copy of the longevity-associated variant. (one copy, evidence: moderate)

### Metabolic
- **ALDH2** (Alcohol flush): Typical alcohol processing. (no copies, evidence: strong)
- **FTO** (Body weight set point): Slightly higher body weight set point at this locus. (one copy, evidence: strong)
- **TCF7L2** (Type 2 diabetes risk): Modestly raised type 2 diabetes risk at this locus. (one copy, evidence: strong)

### Pharmacogenomics
- **SLCO1B1** (Statin transport): Typical statin transport. (no copies, evidence: strong)

## What this file cannot show

- A consumer array assays a tiny, preselected fraction of the genome. A normal result here never rules out other variants in the same gene.
- 23andMe chip content changed between v3, v4, and v5. A missing marker here usually reflects which chip you were tested on rather than a real negative.
```

Paste that into Claude or ChatGPT. Keep the other 640,000 variants to yourself.

APOE is absent from that output on purpose: it needs `--include-sensitive`, because Alzheimer's risk is a thing you decide to learn rather than encounter while reading about cholesterol.

## Why this exists

The obvious way to ask an LLM about your DNA is to upload the raw file. That file has around 640,000 positions in it. It identifies you permanently, it reveals things you were not asking about (carrier status, disease risk, whether your father is your father), and it does the same for your siblings and children, who were not asked.

Almost none of that is needed to have a conversation about how you age.

DNAdistiller reads the file on your machine, pulls out a curated set of variants that bear on longevity and healthspan, and writes a profile you can read in thirty seconds before deciding whether to share it. You choose how much detail goes in.

## What "longevity" means here, and what it excludes

Longevity is not the theme. It is the **entry requirement**, and it is the reason the catalog is twelve markers rather than two hundred.

A marker earns a place only if it plausibly bears on **how long you stay healthy**. In practice that means one of:

- **The things that actually kill people.** Cardiovascular disease and diabetes are most of it, so LPA, TCF7L2, and FTO are in.
- **Common, serious, and fixable.** HFE haemochromatosis is the best marker here by a distance: about 1 in 271 Northern Europeans is a C282Y homozygote, untreated iron overload causes cirrhosis and diabetes and heart failure, and the treatment is donating blood.
- **Decisions made in advance or not at all.** Factor V Leiden and prothrombin G20210A change what a prescriber does about oestrogen contraception, HRT, surgery, and pregnancy. Nothing in a blood panel warns you first.
- **Cognitive ageing.** APOE, which is the largest common-variant effect in the whole catalog.
- **Drug safety over a lifetime.** SLCO1B1, because statins are the drug most people take longest, and a muscle side effect is a reason people stop taking one that was working.
- **Exposures that compound over decades.** ALDH2, where drinking despite the variant raises oesophageal cancer risk substantially.
- **Longevity itself.** FOXO3, the one "longevity gene" besides APOE that survived replication.

That bar is what keeps this a healthspan tool instead of a trait grab-bag. Most consumer DNA reports are grab-bags: they will tell you about earwax type, cilantro aversion, and "warrior vs worrier" personality genotypes, all rendered with the same confidence as a real drug-safety finding. Those are fun. They are not longevity, and they crowd out the handful of results that matter.

**Two markers sit at the edge and are worth naming:**

- **CYP1A2** (caffeine) is in on a healthspan argument, not a "know thyself" one: coffee intake interacts with this genotype on heart-attack risk (Cornelis 2006). It stays because of that link, not because caffeine is interesting.
- **MTHFR** is in on no scientific merit at all. It is graded `weak`, hidden by default, and exists so that someone who has been told online they have an "MTHFR mutation" can read what the evidence actually says. It is a rebuttal, not a finding.

The catalog is small on purpose: twelve markers at v0.1, each checked against dbSNP and a paper, graded, and cited. Twelve solid markers beat two hundred speculative ones, and most tools in this space chose the two hundred.

## Genetics is a prior. Go get the measurement.

The most useful thing this tool does is tell you which blood tests are worth asking for.

DNA is bad at "how am I doing" and genuinely good at "what is worth measuring". A genotype gives you a lifetime prior; a ferritin gives you the actual number, today, for about the price of lunch. So every profile ends with a **Tests worth asking for** section, split three ways:

- **The test beats the genotype.** Most of the catalog. Ferritin over HFE, because iron loading is what harms you and most C282Y homozygotes never load iron. HbA1c over TCF7L2. Homocysteine over MTHFR, which will probably come back normal, and that is the point.
- **Worth having alongside.** It answers a different question than the variant.
- **No test warns you in time.** The small, valuable set: nothing tells you your Factor V Leiden status before the clot, and creatine kinase only rises after the statin has already damaged muscle. These are the markers where the genotype *is* the point.

A tool that reports the prior and never names the measurement has stopped halfway. Most of this catalog defers to a blood test, and there is a test asserting that stays true.

## What honest looks like here

This tool will tell you that most of your results do not matter very much. That is the main feature.

Consumer genomics has a credibility problem, and it comes from tools that take an odds ratio of 1.12 and report it as "elevated risk". So:

- Every marker carries an evidence grade, and markers graded weak are hidden unless you ask for them. MTHFR is in the catalog, graded weak, with the evidence against it written down. It is there because leaving it out looks like an oversight, not because it means much.
- Markers the chip never tested are listed as "not tested" instead of quietly left out. A gap in a risk list reads as an all-clear, and an array cannot support that. It tests preselected positions. It does not sequence.
- The profile ships with a preamble aimed at the LLM, telling it to stay calibrated and not to convert small effects into predictions. Without it, models reliably over-read this data.
- Effect sizes are printed next to the findings, so an OR of 1.08 looks like what it is.

## Install

Requires Python 3.10 or newer.

```sh
uv tool install dnadistiller
```

Or from source:

```sh
git clone https://github.com/Lucaks3/DNAdistiller
cd dnadistiller
uv venv && uv pip install -e ".[dev]"
```

## Getting your DNA file

Download the raw data export from your provider. It is a text file, usually 15-25 MB unzipped.

- 23andMe: Settings, then "23andMe Data", then Download raw data
- AncestryDNA: Settings, then DNA, then Download raw DNA data
- MyHeritage: DNA menu, then Manage DNA kits, then Download

Also supported: FamilyTreeDNA, Living DNA, and single-sample VCF exports.

## Usage

Check that your file parses before anything else:

```sh
dnadistiller check ~/genome.txt
```

This reads no markers and interprets nothing. It tells you which provider it detected, how many variants it found, the reference build, your no-call rate, and how much of the catalog your particular chip can answer. If this step is wrong, everything after it is wrong in ways that are hard to spot.

Build a profile:

```sh
dnadistiller profile ~/genome.txt                          # standard tier, sensible defaults
dnadistiller profile ~/genome.txt --tier minimal           # share less
dnadistiller profile ~/genome.txt -c caffeine -c lipids    # only these topics
dnadistiller profile ~/genome.txt -o profile.md            # write to a file
dnadistiller profile ~/genome.txt -f json                  # for piping elsewhere
```

See what the tool claims before you trust it with anything:

```sh
dnadistiller markers --show-weak
```

## Disclosure tiers

The same finding, at each tier:

| Tier | Output |
|---|---|
| `minimal` | `Intermediate caffeine metabolism.` |
| `standard` | `**CYP1A2** (Caffeine metabolism rate): Intermediate caffeine metabolism. (one copy, evidence: moderate)` |
| `full` | `**CYP1A2** (rs762551): AC (one copy)` plus effect size, ancestry caveat, and citations |

`standard` is the default. `full` includes rsIDs and genotypes, which is what you want if you plan to check a finding against the literature, and it is the most identifying option.

`minimal` drops the gene name and the category too. That costs you something real: the model can no longer look the variant up. It is the right choice when you want to talk about the finding rather than the variant.

Markers whose results are hard to un-know are excluded from every tier unless you pass `--include-sensitive`. Right now that means APOE, where the ε4 allele carries most of the common genetic signal for late-onset Alzheimer's. Genetic counsellors treat that disclosure as a decision a person makes on purpose. So does this tool.

## What "safer" does and does not mean

Sharing twelve markers instead of 640,000 is data minimisation. It is not anonymisation, and the difference matters.

Genotypes at a few dozen independent common variants are enough to single a person out. A `full` profile is therefore still identifying, and a `standard` one is still information about your body that you cannot take back once sent. What the tiers buy you is a choice about how much you hand over, and a ceiling on the damage when the recipient turns out to be less careful than you assumed.

Before pasting anything into a chatbot, know that provider retention and training policies differ by tier and change over time. Consumer chat products have historically defaulted to training on your conversations. `docs/privacy.md` goes into this properly.

## What this tool will not do

It will not diagnose anything, and it is not a medical device. Nothing here is medical advice.

It will not tell you your ancestry, find relatives, or compute a polygenic risk score. Those need the whole file, which is the thing this tool exists to avoid sending anywhere.

It will not interpret rare or clinically actionable variants like BRCA. Consumer arrays have false-positive rates at rare positions high enough that a single result is not trustworthy, and the failure mode is somebody making a decision on a bad call. That needs a clinical-grade test and a genetic counsellor, not a CLI.

It will not upload your data, because it has no code that can. Read `src/dnadistiller/` and check: roughly 2,500 lines across eleven files, a good share of which is comments explaining why, and three dependencies (`typer`, `rich`, `pyyaml`). CI fails the build if anyone adds an import that could open a socket. That is the point of keeping it small. A promise you can verify in an afternoon is worth more than one you have to take on faith.

## Contributing

Markers live in `data/markers/*.yaml` as data, not code. Adding one is a small YAML file with an rsID, the genotypes, and citations. You do not need to read Python to review whether a claim is true, which is the point of the split.

`CONTRIBUTING.md` has the marker schema and what a claim needs before it gets graded strong or moderate.

Never commit a real DNA file, including your own. `tests/fixtures/README.md` explains why the rule has no exceptions.

## Licence

Apache 2.0. See `LICENSE`.

## Disclaimer

For education and research. Not a medical device, not a diagnostic, not a substitute for a clinician. Consumer genotyping arrays are not clinically validated. Do not make medical decisions based on this output. If something here concerns you, take it to a doctor or a genetic counsellor, who will start by ordering a test that is actually validated.
