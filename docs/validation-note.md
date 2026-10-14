# Validating a diet–vaginal-microbiome analysis toolkit on public cohorts

**Jean Yu · October 2026 · toolkit:** [github.com/jyu9675/vmbdiet](https://github.com/jyu9675/vmbdiet)

## Aim

Several studies report that diet shapes the vaginal microbiome, but each uses its own ad-hoc pipeline, which
makes the signal hard to compare or build on. I packaged the analysis the field shares — community state type
(CST) classification, then covariate-adjusted association of diet with CST / BV / *Lactobacillus* dominance —
into one tested, reusable tool (`vmbdiet`), and asked a concrete question: **does it recover a known published
result when run, from raw data, on a cohort it was not tuned to?**

## Approach

Two fully public cohorts that pair vaginal 16S with dietary intake, each **reprocessed from the raw SRA reads**
(no reuse of the authors' processed tables), with one open-source pipeline: cutadapt → vsearch pair-merge →
UNOISE3 ASVs → RDP genus + curated NCBI-RefSeq 16S for species-level *Lactobacillus* → CST by the dominant-taxon
rule (VALENCIA centroids also shipped). Associations — per-SD adjusted logistic regression with BH-FDR, Spearman
nutrient↔taxon correlation, and PERMANOVA — are computed by `vmbdiet.associate`. Everything runs on a laptop.

## Results

**Cohort 1 — Italian "StressDiet" (PRJNA1188525, N = 113, non-pregnant).** The reconstructed communities are
textbook (*L. crispatus* 41%, *L. iners* 22%, *G. vaginalis* 6.3%; 31% CST-IV). **Both published diet signals
replicate independently:**
- Animal-protein fraction → CST-IV: energy-adjusted **OR 1.66 per SD (1.09–2.53), p = 0.019**.
- Alcohol → *Gardnerella*: Spearman **ρ +0.30, p = 0.001**; *Leptotrichia* ρ +0.36 survives FDR (q = 0.012).

Energy adjustment was decisive: on raw intake, fiber and vegetable protein looked protective, but that was a
confound of lower total energy in CST-IV women; nutrient-density modeling removed the artifact and left the
animal-protein effect — i.e. the tool handled a standard nutritional-epidemiology pitfall correctly.

**Cohort 2 — MELODY (PRJNA915128, N = 47, pregnancy ± IBD), HEI-2015.** A concordant-direction **null**: no diet
variable reaches significance (animal protein OR 1.30, p = 0.44; HEI-2015 OR 1.40, p = 0.33). This matches the
source study's own conclusion of no diet–microbiota association, and is expected — the cohort is underpowered
(12 CST-IV events), pregnancy drives *Lactobacillus* dominance regardless of diet, and alcohol intake is ≈ 0 so
that axis cannot be tested. The one comparable signal (animal protein) still points the same way as Cohort 1.

## Interpretation

A positive replication in the powered cohort and an honest null in the underpowered one is the behavior a
trustworthy tool should show — it does not manufacture significance where there is no exposure contrast or
sample size. The animal-protein→CST-IV direction is consistent across both cohorts.

**Limits:** both cohorts are modest and cross-sectional; the animal-protein association is nominally significant
but does not clear FDR across all nutrients; species resolution is limited to the dominant vaginal taxa.

## What this enables

The pipeline and toolkit are now validated end-to-end and reusable: the same code can profile any cohort that
pairs a 16S table with diet, which makes a **prospective or larger-cohort diet→CST study**, or a multi-cohort
meta-analysis, a low-risk next step with the analysis already in hand. It also complements the lab's existing
strain-transmission / engraftment work — diet is a plausible modifier of which *Lactobacillus* strains establish
and persist. The immediate open item is that the Hawaii multi-ethnic HEI-2015 cohort (PMC11479099) has not yet
deposited an accession; it would be the ideal third HEI-2015 replication once public.

*Reproducible scripts, derived tables, and the full write-up:
[`example/`](../example/) and [`docs/italy-stressdiet-findings.md`](italy-stressdiet-findings.md).*
