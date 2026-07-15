---
name: Propose a marker
about: Suggest a new locus for the catalog
title: "[marker] GENE rsID"
labels: marker
---

<!--
Never paste your own genotype or a raw DNA file into an issue. Describe the
marker, not your result.
-->

## The marker

- **Gene:**
- **rsID:**
- **Category:** (lipids / metabolic / methylation / caffeine / inflammation / neuro / longevity / fitness / sleep / pharmacogenomics)

## The evidence

- **Primary citation (PMID or DOI):**
- **Study size and design:**
- **Has it replicated in independent cohorts?**
- **Effect size:** (the actual number: OR, beta, or HR. "Increased risk" is not an effect size.)
- **Populations studied:** (and whether it is known to transfer)

## Proposed grade

- [ ] `strong`: large replicated GWAS or meta-analysis, or pharmacogenomics with dosing guidelines behind it
- [ ] `moderate`: replicated, modest effect
- [ ] `weak`: small, inconsistent, or unreplicated

## Alleles

- **Effect allele:**
- **Strand:** Have you checked the plus-strand alleles on [dbSNP](https://www.ncbi.nlm.nih.gov/snp/)? This is the most common way a marker ends up silently wrong.
- **Present on which chips?** (23andMe v3/v4/v5, AncestryDNA, unknown)

## Why it belongs here

What does knowing this change for someone? "It is interesting" is a fair answer, but say so rather than implying it is actionable.

## Sensitivity

- [ ] This finding is hard to un-know (like APOE and Alzheimer's) and should be gated behind `--include-sensitive`
