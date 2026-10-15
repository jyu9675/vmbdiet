# MaAsLin 3 bridge

`vmbdiet.maaslin` adds the repeated-measures association model the gut–vagina-axis / FRESH design
requires. The existing `vmbdiet.associate` module tests one sample at a time (fixed-effects logistic
OR, Spearman) — correct for a single vaginal swab per person, but wrong when each **participant**
contributes **rectal + vaginal** samples across **multiple timepoints**, because those samples are
correlated. The field-standard answer is [MaAsLin 3](https://bioconductor.statistik.tu-dortmund.de/packages/3.22/bioc/vignettes/maaslin3/inst/doc/maaslin3_tutorial.html),
which fits a per-feature **abundance** model and a separate **prevalence** model and absorbs the
within-host correlation with a random intercept `(1 | participant)`.

## Why two pieces

MaAsLin 3 is R/Bioconductor only. This module therefore provides:

1. **`maaslin_py(...)`** — a faithful Python re-implementation that runs with only `statsmodels`:
   - *abundance*: `log(relative abundance) ~ predictors + covariates [+ body_site]` with a random
     intercept for the grouping column, fit on the **present** samples (mixed LM; OLS fallback when the
     random-effect variance is unidentifiable).
   - *prevalence*: `present(0/1) ~ same fixed effects` as a logistic **GEE** with exchangeable
     within-group correlation (clustering-aware).
   - BH-FDR within each model; per-SD effects when `standardize=True`; per-feature failures captured in
     a `note` column, never raised.
   This is what lets the whole pipeline be developed, tested, and demoed on a laptop **before** the real
   HUMAnN 3 / species tables or controlled-access cohorts are available.

2. **`prepare_inputs(...)` + `write_maaslin3_script(...)`** — emit MaAsLin 3-format inputs and a
   ready-to-run R script with the correct formula, for the lab cluster. Built from the **same**
   arguments as the Python run, so the two fits are directly comparable.

## Use

```python
from vmbdiet import maaslin_py
res = maaslin_py(feature_table, metadata,
                 predictors=["animal_protein", "fiber"],
                 covariates=["age", "race", "antibiotics"],
                 group="participant", body_site="body_site")
```

or from the CLI (`vmbdiet maaslin --features ... --meta ... --predictors ... --group participant
--body-site body_site`). The feature table is species relative abundances (Aims 1–2) **or** HUMAnN 3
pathway / eggNOG gene-family abundances (Aim 3) — identical call. See
[`example/maaslin_demo`](../example/maaslin_demo) for an end-to-end synthetic run.

## Scope / what it is not

- The Python engine mirrors MaAsLin 3's *structure* (two models, random intercept, TSS+log, BH); it is
  **not** a byte-for-byte reproduction of MaAsLin 3's internals (e.g. its median-comparison abundance
  null, augmentation, or variance weighting). For reported results, run the generated `run_maaslin3.R`;
  use the Python engine for development, CI, and sanity-checking the R output.
- **Mediation is out of scope here.** The axis hypothesis — diet → gut β-glucuronidase/estrobolome →
  vaginal glycogen → *L. crispatus* — is a *mediation* chain; MaAsLin 3 (and this bridge) test only the
  marginal links (diet↔function, diet↔CST). A mediation model is the natural next module.
