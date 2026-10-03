# vmbdiet — design note

## The question

Vaginal community composition is a strong correlate of reproductive and pregnancy health: a
*Lactobacillus*-dominant community (CST I/II/III/V) is protective; a diverse, anaerobe-rich community
(CST IV) is bacterial vaginosis and carries elevated risk of preterm birth, STI and HIV acquisition.
Composition is usually treated as given. But several cohorts find it tracks a **modifiable exposure —
diet**. `vmbdiet` standardizes the test of that association so it runs identically on any cohort that
pairs a 16S table with dietary intake.

## The analysis, as the literature does it

Every diet–vaginal-microbiome study runs the same two-step analysis; they differ mainly in cohort and
in which nutrients they test. `vmbdiet` implements the shared core:

1. **Classify the community.** Assign each sample a CST (dominant *Lactobacillus* species, or diverse/BV).
   `cst.assign_cst` uses the transparent dominant-taxon rule (dominant taxon ≥ 50%); `assign_cst_valencia`
   reproduces VALENCIA's nearest-centroid (Yue-Clayton) assignment when the reference centroids are supplied.
2. **Associate diet with the community**, adjusting for confounders:
   - **Adjusted logistic regression** (`logistic_assoc` / `logistic_scan`): `P(outcome) ~ nutrient + covariates`,
     where outcome ∈ {BV, CST-IV, *Lactobacillus*-depleted}. Reported as a **per-SD odds ratio with 95% CI**,
     with Benjamini-Hochberg FDR across nutrients. This is the primary model in Neggers 2007 and the mixed/
     adjusted models in the prenatal and HEI studies. (A random intercept for repeated visits is the natural
     extension for longitudinal cohorts.)
   - **Nutrient ↔ taxon correlation** (`nutrient_taxon_corr`): Spearman/Pearson of a nutrient against each
     taxon's abundance (the "carbohydrate ↔ *L. crispatus*" style result), FDR-corrected.
   - **Mann-Whitney** (`mannwhitney_by_group`): nutrient intake, *Lactobacillus*-dominant vs depleted.
   - **PERMANOVA** (`permanova` on a Bray-Curtis matrix from `diversity`): does overall composition differ
     across a diet grouping — the beta-diversity test (PCoA via `diversity.pcoa`).

## Design choices

- **Per-SD odds ratios.** Standardizing continuous nutrients makes ORs comparable across nutrients on
  different scales (g/day vs µg/day) — you can rank effects directly. Raw-unit ORs are available (`standardize=False`).
- **Covariate adjustment is first-class.** Every study adjusts for age, race/ethnicity, BMI, and often
  alcohol, contraception, income/education. Categorical covariates are auto-wrapped as factors.
- **FDR, not raw p.** Nutrient panels are wide; BH-FDR is applied across the scanned nutrients and across taxa.
- **Dependency-light CST default.** The dominant-taxon classifier needs only a species-named abundance
  table — no reference file — so the tool runs anywhere; VALENCIA centroids are an opt-in upgrade for exact
  sub-CST alignment with the literature.

## Validation

- **Synthetic (done, `tests/`)** — data generated with fat injected as a CST-IV risk factor and
  carbohydrate/folate as protective; the tools recover exactly that (fat OR > 1, carb/folate OR < 1, all
  FDR-significant; PERMANOVA p ≈ 0.001). Proves the statistics and the CST logic.
- **Real public data (next)** — the Italian cohort (PRJNA1188525: 16S + nutrient FFQ) is the first target;
  MELODY (PRJNA915128, HEI-2015) the HEI cross-check. See [DATASETS.md](DATASETS.md).

## Limits (v0.1)

- Dominant-taxon CST is a proxy for VALENCIA's sub-CSTs; use `assign_cst_valencia` with the published
  centroids for exact sub-CST labels.
- Logistic models are cross-sectional; longitudinal cohorts (repeated visits) want a mixed-effects model
  with a per-subject random intercept (planned).
- Association is not causation — these are observational diet exposures with the usual confounding caveats,
  which is why covariate adjustment is built in rather than optional.
