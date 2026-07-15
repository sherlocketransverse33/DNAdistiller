# Privacy: what this tool does and does not buy you

The pitch is "share 40 markers instead of 640,000". This page is the honest version of that claim, including the parts that undercut it.

## Minimisation, not anonymisation

**A profile from this tool is not anonymous.** Not at any tier. If you take one thing from this page, take that.

Genotypes at a few dozen independent common variants are enough to distinguish one person from everyone alive. The number has been in the literature since Lin, Owen and Altman put it at roughly 30 to 80 SNPs in 2004 ([Science 305:183](https://doi.org/10.1126/science.1095019)), and a `full` profile is in that territory by construction. Sharing one is closer to handing over a fingerprint than to sharing a statistic.

What the tool actually does is reduce how much you disclose and give you a say in it. That is worth real money in a way "anonymous" would not be:

- **Less surface.** A raw file carries ancestry, relatedness, carrier status for conditions you never asked about, and every risk marker on the chip. A longevity profile carries the markers you asked for.
- **Fewer surprises.** You cannot accidentally disclose your APOE status if APOE was never in the file you sent.
- **A ceiling on the damage.** If the recipient turns out to be less careful than you assumed, what they have is bounded by what you chose.

That is data minimisation. It is a real privacy improvement and it is the standard defensible claim. Anyone marketing a subset of your genome as "anonymised" is either confused or selling something.

## What each tier actually discloses

| Tier | Contains | Realistic risk |
|---|---|---|
| `minimal` | Interpretations only. No rsIDs, genotypes, gene names, or categories. | Lowest. A recipient cannot look the variants up or match you against a reference database from this alone. It still describes your body. |
| `standard` | Gene names and copy counts. | Moderate. Gene plus copy count is recoverable to something close to a genotype for anyone who knows the biology. Not lookup-ready, but not opaque. |
| `full` | rsIDs and genotypes. | Highest, and identifying. Treat it like a medical record. |

`standard` is the default because it is useful and bounded. `full` exists because checking a finding against the literature needs the rsID, and that is a legitimate thing to want.

### On `minimal`

Stripping rsIDs and genotypes genuinely helps: "slow caffeine metabolism" is not a database key, and you cannot run a match against a reference panel with it.

It is not magic. Enough interpretations still narrow a person down, and a rare finding narrows them a lot. `minimal` is the least you can share and still have a conversation worth having. It is not nothing.

## Your relatives did not agree to this

Your genome is substantially your parents', siblings', and children's genome. Sharing yours discloses theirs, and they were not asked.

This is why the project treats "it's my data, my choice" as incomplete rather than wrong. It is your choice. It is not only your data.

## Before you paste anything into a chatbot

Sending a profile to an LLM sends it to a company.

Consumer chat products have historically defaulted to training on your conversations, with an opt-out you have to find. Business, enterprise, and API tiers usually do not train on your data by default. These policies differ per product and change, so check yours today rather than trusting what was true last year, and check whether retention has a floor even after you delete a chat.

If that gives you pause, the honest options are to use `--tier minimal`, to use a tier of the product that does not train on your input, or to run a local model. All three are reasonable. So is deciding this is not worth doing.

## The thing that actually creates risk is a diagnosis, not a file

Even Florida's ban, the strongest genetic-insurance law in the US, has a trapdoor worth understanding, because it generalises.

Section 627.4301(2)(a) protects you *in the absence of a diagnosis*. Then (2)(d) expressly permits life and long-term-care insurers to consider *"a medical diagnosis... even if [it] was made based on the results of a genetic test."*

Read those together and the protection inverts. It holds while your genotype is just a genotype. It evaporates the moment that genotype becomes a diagnosis written into your medical record. The strongest genetic-privacy law in the country stops protecting you at exactly the point a doctor writes something down.

Which sharpens the whole point of this tool:

> The risk is not the file on your laptop. It is what you tell an insurer, and what your doctor writes in your chart.

A local tool whose output never enters a medical record is a genuinely different risk profile from a clinical test, and that is the strongest honest thing to say for local-first. It is also not a reason to avoid a doctor. If a result here matters, it is worth a real test and a real conversation, and that conversation is worth more than the insurance consequence in almost every case. Just make the trade knowingly rather than by accident.

## Why local-first is not paranoia

23andMe filed for bankruptcy in 2025, and the question of who would end up owning fourteen million people's genomes became a live one. Nobody who spat in that tube in 2013 was told that was a thing they were agreeing to.

That is the structural problem with custody: policies are promises by an entity that may not outlive your data. Genetic data is permanent, and consent given to one company is worth whatever a bankruptcy court says it is worth.

This tool has no servers, no accounts, and no network code, which is why. Not because a server would be badly run, but because the safest custodian of your genome is the machine it is already on.

## The law, roughly

Not legal advice, and the two hardest questions below have no published authority behind them at all. Jurisdictions differ and this changes.

- **GDPR** treats genetic data as an Article 9 special category, with a higher bar than ordinary personal data. Software that runs locally and transmits nothing largely sidesteps this: no controller, no processing, no transfer. Pasting a profile into a chatbot creates a transfer, and you are the one making it.
- **GINA** (US) bars genetic discrimination in health insurance and employment. It does **not** cover life insurance, disability insurance, or long-term care insurance. Those are the ones where a disclosed APOE ε4 status can actually cost you money, and they are exactly the ones GINA leaves out.

  Do not lean on state law to fill that gap without reading the statute. The protections are much thinner than the trackers suggest: one widely circulated chart claims roughly 46 states have life-insurance genetic protections, and North Dakota and Texas turn out to have nothing at all. Florida is the only true flat ban for life insurance. Maryland joined in October 2025 as an actuarial hybrid rather than a ban. Nevada, Wisconsin, and Kansas get counted as "genetic privacy states" while exempting the lines that matter. The error rate in secondary sources here is remarkable. Check the statute.
- **HIPAA** almost certainly does not apply. It binds covered entities like providers and insurers, not a CLI on your laptop or a chatbot you paste into.

### Medical device rules bite harder than you would guess

Three intuitions that feel like defences and are not:

- **"It runs locally."** MDR Recital 19 says qualification is *independent of the software's location*. Local execution is a privacy answer, not a regulatory one.
- **"It's free."** MDR Article 2(27) covers supply *"whether in return for payment or free of charge"*. Price is expressly irrelevant.
- **"It examines no specimen, so IVDR cannot apply."** This is the one that surprised us. MDCG 2019-11 tests the *data source*, not specimen contact, and its Annex I lists as an example of IVD software: *"MDSW that integrates genotype of multiple genes to predict risk a disease or medical condition developing or recurring."* That is a description of this tool. Generating new information from already-available genotype data is enough.

A fourth one is worth killing carefully, because it is the intuition this file previously got wrong: **"there's a disclaimer on it."** Under 21 CFR 801.4 intended use includes the product's design and circumstances, not only what you say about it. Marketing can create device status, and FDA's 2013 warning letter to 23andMe rested entirely on marketing copy. Disclaiming cannot remove it. A "not medical advice" banner does not fix a tool whose architecture exists to surface disease risk.

What actually moves the needle is narrower and less comfortable:

> Name a disease and every regime engages at once. Name none and most never start.

That is a product decision, not a legal one, and it converges from several directions: 21 CFR 801.4, FDA's general-wellness decision tree (which a genotype tool fails twice, since any disease reference fails the first test and a genotype is not a lifestyle choice), Germany's GenDG § 3 Nr. 8 with its "disease arising only in the future", and MDR Article 2(1), which added "prediction" and "prognosis" that the old directive never had.

**So be clear-eyed about where this project sits.** It names Alzheimer's, type 2 diabetes, and coronary disease, because saying "APOE relates to a condition we decline to name" would be worse for you and no safer for us. The wellness carve-out in MDR Recital 19 is therefore not a defence this tool can honestly claim. What it has instead is a set of choices that keep the exposure as small as it can be while still being useful: no diagnostic claims, no dosing advice, no BRCA-class variants, a preamble that tells the model to defer to a clinician, and a catalog that describes rather than prescribes.

Those are real, and they are not a safe harbour. Anyone forking this to make it prescriptive should get a notified-body view first.

### The part that lands on you, not on the project

In several European countries the binding constraint is national law, and it is aimed at the person running the test rather than the person who wrote the software.

- **Germany.** The Gendiagnostikgesetz reserves predictive genetic examinations to qualified physicians (§ 7), with counselling before and after (§ 10). Whether that reaches a tool that performs no laboratory analysis is genuinely unsettled and there is no published authority either way. If it does, the person performing the examination is the user, at their own keyboard, and the exposure is an administrative fine of up to €50,000.
- **France.** Code civil article 16-10 permits examination of genetic characteristics only for medical or research purposes. Code pénal article 226-28-1 then penalises *the person who solicits examination of their own genetic characteristics* outside those conditions, at €3,750. France fines the consumer, not the provider. Reinterpreting a file you already hold arguably solicits nothing, so this is likely weaker than Germany, but nobody has tested it.

We flag this rather than hiding it because the asymmetry is the uncomfortable part: the risk, such as it is, sits with the user and not with the person who published the code. That is worth knowing before you run it.

## Why the tool refuses to interpret clinically actionable variants

dnadistiller will not report BRCA or similar.

Consumer arrays are not clinically validated, and at rare positions their false-positive rate is high. A 2018 study in *Genetics in Medicine* ([Tandy-Connor et al.](https://doi.org/10.1038/gim.2018.38)) found that around 40% of variants flagged in consumer raw data were false positives on confirmatory testing.

Forty percent is not a caveat, it is a coin flip. The failure mode is somebody making an irreversible decision on a bad call from a spit tube. That needs a clinical-grade test and a genetic counsellor, and no amount of disclaimer text in a CLI substitutes for either.

## Threat model

**Defends against:** disclosing far more than you meant to, to a party you chose, because the only available unit of sharing was "everything".

**Does not defend against:** anyone who already controls your machine. If they can read your home directory they can read the DNA file itself and have no need of us.

**Cannot defend against:** what the recipient does next. Once a profile is sent, it is theirs. The tiers decide what "it" was.
