# Architecture and the decisions behind it

Written down because each of these looks like an oversight until you know the reason, and someone will otherwise "fix" them.

## Lookups key on rsID, never on position

A marker is found by its rsID. Chromosome and position are parsed, reported by `dnadistiller check`, and then never used to find anything.

Positions move between genome assemblies. rs429358 sits at 19:50103781 in NCBI36 and 19:45411941 in GRCh37, about 4.7 Mb apart. A tool that looks up markers by coordinate has to know the build, get the remapping right, and keep doing so as assemblies change. When it gets that wrong it does not fail: it reads a different locus and reports a confident answer about the wrong variant.

Keying on rsID makes a build mismatch impossible rather than merely unlikely. The cost is that a variant without an rsID cannot be catalogued, which has not been a real constraint.

This is also why build 36 files from 2011 still work. Build is provenance we report, not a coordinate system we have to reconcile.

## No pandas

The original prototype imported pandas to call `value_counts()` on a 600,000-row file.

The actual work is looking up around 40 known rsIDs. That is a dict built from a line iterator. Plain Python does it faster, in constant memory, and without a 50 MB dependency.

This matters more than it would elsewhere. The tool's claim is that you can audit it before trusting it with your genome. Three dependencies and 900 lines is a claim someone can check in an afternoon. Numpy and pandas underneath is not.

## We do not use the `snps` library

[`snps`](https://github.com/apriha/snps) is a good library. It is BSD-3-Clause, actively maintained, reads about twenty formats to our six, and has already encountered every parsing quirk we have. Adopting it would be the obvious call, and we did not.

Two reasons:

**Its main value does not apply to us.** The hardest thing `snps` solves is remapping positions between builds. We do not remap positions, because we do not use positions. That entire capability is dead weight here.

**Its dependencies break the one promise the project makes.** `snps` pulls in `pooch`, a data downloader, and does build remapping by calling the Ensembl API. Both are network paths. "There is no code here that can send your genome anywhere" stops being verifiable the moment the dependency tree can open a socket, and that sentence is the whole product. A user cannot audit a promise that rests on a transitive dependency never being invoked.

So we own our parsers, and we own their bugs. The mitigation is that the design makes parser bugs loud instead of silent: an unrecognised call reports `UNKNOWN_GENOTYPE`, an absent rsID reports `NOT_ON_CHIP`, and neither is ever rendered as "normal".

If you need build remapping, exotic formats, or ancestry work, use `snps`. It is the better tool for that job. This one has a narrower job.

## The full genome never touches disk

`Genome` holds every call and lives in memory only. Nothing serialises it. The only object designed to be written out is `Profile`, which is redacted by construction.

The prototype this replaced uploaded the file to a Flask server, saved it to a fixed path, and deleted it only on the success path, so a parse error left genetic data on disk indefinitely. That is the failure this rules out structurally rather than carefully.

## Redaction happens in exactly one function

Every output format calls `profile.redact`. A renderer receives a dict of permitted fields and physically cannot disclose one it was not handed.

The alternative, where each renderer reads `Finding` and remembers what its tier allows, fails the first time someone adds a format. During development the JSON renderer disclosed a field the Markdown renderer withheld, at the same tier. The test that caught it asserts on rendered output rather than on internal state, because what matters is what ends up in the file someone pastes.

## Filters are opt-in

`select_markers` starts from nothing and adds what was asked for. Sensitive markers, weak evidence, and unrequested categories are never in the set to begin with.

A bug in an opt-out filter discloses more than the user wanted. The same bug in an opt-in filter discloses less. Only one of those is recoverable, and the user does not find out about the other one.

## The preamble is part of the product

`prompts.PREAMBLE` is aimed at the LLM, not the user, and it is emitted by default.

Pasting a bare list of risk alleles into a chatbot reliably produces a confident narrative built on odds ratios of 1.1. The preamble measurably changes that answer. It is not a disclaimer protecting us; it is the difference between the tool being useful and being actively misleading, so it is not optional and it redacts along with everything else.

## Strand is resolved at interpretation, not at parse

`Genotype.alleles` holds what the provider wrote. Normalisation happens in `interpret._match`, which tries a direct match and falls back to the complement.

The correct strand depends on which marker you are asking about, so there is no strand to normalise *to* at parse time. Doing it later means the catalog declares what it expects and the check is local to the marker.

The fallback is safe because it is asymmetric: a genotype that already matches the catalog is by definition already on the catalog's strand, so the complement is only ever tried after a miss. It cannot work at A/T or C/G sites, where the complement of a valid genotype is another valid genotype. There the tool declines and reports the marker unread, because a coin flip between two readings is worse than no reading.

## Every failure to read a marker is reported

`NOT_ON_CHIP`, `NO_CALL`, and `UNKNOWN_GENOTYPE` are distinct, and all three appear in output under "Not tested".

A list of genetic findings with silent gaps reads as an all-clear on everything absent. A genotyping array tests preselected positions and says nothing whatsoever about the rest, so that reading is unsupportable and the tool must not invite it.
