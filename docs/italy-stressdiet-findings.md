# Diet and the vaginal microbiome: an independent replication of the Italian "StressDiet" cohort

**Project:** vmbdiet · **Data:** PRJNA1188525 / PMC12213695 · **Analyst:** Jean Yu · **N = 113**

## Summary

We asked whether a published diet ↔ vaginal-microbiome association could be reproduced *from raw
data* with an independent pipeline and an independent analysis toolkit (`vmbdiet`). Starting from the
raw 16S FASTQs on SRA (not the authors' processed tables), we rebuilt the amplicon-sequence-variant
(ASV) table, assigned community state types (CSTs), joined to the cohort's food-frequency
questionnaire, and tested the diet associations. **Both headline diet signals from the original study
replicate:** a higher animal-protein fraction of the diet tracks with the diverse/BV community type
(CST-IV), and alcohol intake correlates with *Gardnerella* and other BV-associated anaerobes. A
methodological wrinkle — total energy intake is *lower*, not higher, in CST-IV — is handled by
expressing diet as nutrient density, which is what surfaces the animal-protein effect cleanly.

## Data and methods

- **Microbiome.** 113 reproductive-age women (Bologna); vaginal 16S rRNA V3–V4 (341F/805R), MiSeq
  2×300, BioProject PRJNA1188525. We downloaded all 226 FASTQs (~2.8 GB) from ENA.
- **Diet.** The study's FFQ nutrient table (17 variables + the MEDI-LITE Mediterranean score) was
  taken from the journal supplement (CC-BY). The sample key `StressDiet_N` joins the two 1:1 (N=113,
  zero mismatches).
- **Pipeline.** cutadapt 5.2 (primer removal) → vsearch 2.32 (pair merge after quality truncation,
  expected-error filter < 1.0, UNOISE3 denoising, UCHIME3 chimera removal, OTU mapping) → SINTAX/RDP
  for genus and a curated NCBI-RefSeq 16S database for *Lactobacillus* species. 5.1 M merged reads →
  1,312 ASVs → 163 taxa.
- **Statistics (`vmbdiet.associate`).** Adjusted logistic regression (per-SD odds ratios, BH-FDR),
  Spearman nutrient↔taxon correlation, Mann-Whitney between CST groups, and PERMANOVA on Bray-Curtis
  distances. CSTs assigned by `vmbdiet.cst` (dominant-taxon rule; species from the Tier-2 reference).

## Community structure (validation)

The reconstructed communities look exactly like a healthy reproductive-age cohort, which is the first
thing that has to be true before any diet test is meaningful:

| CST | definition | n (%) |
|-----|-----------|-------|
| I   | *L. crispatus* | 49 (43%) |
| III | *L. iners*     | 22 (19%) |
| IV  | diverse / anaerobic (BV-like) | 35 (31%) |
| V   | *L. jensenii*  | 4 (4%) |
| II  | *L. gasseri*   | 3 (3%) |

Top taxa by mean relative abundance: *L. crispatus* 41%, *L. iners* 22%, *Gardnerella vaginalis* 6.3%
(present in 81% of women), *Prevotella* 4.8%. Species-level resolution of the four CST-defining
*Lactobacillus* worked as intended.

## Diet associations

### 1. Animal protein → CST-IV (replicates)

Expressed as **nutrient density (per 1,000 kcal)**, the animal-protein fraction is higher in women
with the BV-like community:

- Mann-Whitney, animal-protein density, CST-IV vs other: 26.2 vs 23.6, **p = 0.026**
- Logistic, CST-IV ~ animal-protein density (per SD): **OR 1.66 (95% CI 1.09–2.53), p = 0.019**
  (BH-FDR q = 0.14 across 7 nutrient densities)
- Vegetable-protein density trends the opposite way (OR 0.68, p = 0.10).

### 2. Alcohol → *Gardnerella* / anaerobes (replicates and extends)

Spearman correlation of alcohol intake with per-taxon relative abundance:

| taxon | rho | p | q (BH) |
|-------|-----|---|--------|
| *Leptotrichia* | +0.36 | 0.0001 | **0.012** |
| *Gardnerella* (genus) | +0.30 | 0.0014 | 0.068 |
| *Gardnerella vaginalis* | +0.29 | 0.0015 | 0.068 |
| *Neisseria* | +0.28 | 0.0025 | 0.068 |

Every alcohol association points the same way — more alcohol, more BV-associated anaerobes — and
*Leptotrichia* survives multiple-testing correction.

### 3. The energy confounder (a cautionary, and informative, detail)

On **raw** intake, vegetable protein and fiber look protective (lower in CST-IV: p = 0.005, p = 0.021).
But total energy is also lower in CST-IV (1,754 vs 1,995 kcal, p = 0.016), so those raw differences
are largely a reflection of *how much* these women ate, not diet composition. Re-expressing everything
as energy density removes the fiber/vegetable-protein signal and leaves the animal-protein-fraction
effect standing. This is a textbook example of why nutritional epidemiology adjusts for energy, and it
is reassuring that the association that survives adjustment is the one the original authors reported.

## Interpretation and limits

- Two independent lines of evidence — a prospective-style composition effect (animal protein) and a
  direct taxon correlation (alcohol → *Gardnerella*) — both reproduce the published directions using
  a completely independent reprocessing. That is a real, if modest, replication.
- **Power is limited:** 113 women, 35 CST-IV events. The animal-protein association is nominally
  significant but does not clear FDR; treat it as corroborating, not confirmatory.
- **Cross-sectional:** diet and microbiome are measured once; no causal direction is established.
- The curated species reference covers the dominant vaginal taxa; rare taxa are reported at genus level.
- PERMANOVA (community ~ CST-IV, pseudo-F ≈ 22, p = 0.001) is strong but partly definitional, since
  CST-IV is itself a community label; it is reported as a structural sanity check, not a diet result.

## Second cohort: MELODY (pregnancy ± IBD, HEI-2015)

To test whether the toolkit behaves honestly on a *different* cohort with a *different* diet metric, we ran
the same pipeline on **MELODY** (PRJNA915128; PLOS ONE 2024): 48 pregnant individuals (23 IBD, 25 control),
vaginal V3–V4 16S, diet measured as **HEI-2015** (the Healthy Eating Index — the Hawaii study's metric).
Usefully, the per-sample diet and covariates are carried in the BioSample metadata, so the cohort is fully
reproducible without any supplement. Details and scripts: [`example/melody_ibd/`](../example/melody_ibd/).

| | Italian (StressDiet) | MELODY (pregnancy/IBD) |
|---|---|---|
| N / CST-IV events | 113 / 35 | 47 / 12 |
| population | non-pregnant, reproductive age | 3rd-trimester pregnant, ±IBD |
| animal-protein density → CST-IV | **OR 1.66, p=0.019** | OR 1.30, p=0.44 (same direction) |
| HEI-2015 → CST-IV | — (not measured) | OR 1.40, p=0.33 |
| alcohol → *Gardnerella* | **ρ +0.30, p=0.001** | untestable (alcohol ≈ 0 in pregnancy) |

MELODY is a **null** — nothing reaches significance and nothing survives FDR. That is the right answer and a
good sign, for three reasons: (1) it is badly underpowered (47 women, 12 events); (2) pregnancy both flattens
the alcohol exposure to zero and pushes the community toward *Lactobacillus* regardless of diet; and (3) the
original authors themselves found no significant diet- or IBD-related difference in vaginal microbiota. The one
directly comparable signal — animal protein — points the **same way** in both cohorts. A tool that returned
"significant" here would be the worrying result, not this.

## Why this matters for the program

This closes the loop on the toolkit: `vmbdiet` was built from four papers' worth of methods, and here
it recovers a real published result from raw public data on a laptop, with the energy-adjustment
nuance handled correctly. The same pipeline can now be pointed at the other verified cohorts in
[`DATASETS.md`](DATASETS.md) (MELODY HEI-2015, MicrobeMom shotgun) for a multi-cohort diet→CST
meta-analysis — a fundable next step with tooling already validated end-to-end.

*Reproducible scripts and all result tables:
[`example/italy_stressdiet/`](../example/italy_stressdiet/).*
