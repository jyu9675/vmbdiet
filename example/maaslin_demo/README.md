# MaAsLin 3 bridge — repeated-measures demo (gut–vagina axis)

A synthetic, laptop-runnable demo of `vmbdiet maaslin`: the repeated-measures association engine the
FRESH / gut–vagina-axis design needs, where each participant contributes **rectal + vaginal** samples
across **two timepoints** (so samples from one host are not independent).

## Run it

```bash
python example/maaslin_demo/make_demo.py          # writes data/features.tsv + data/metadata.tsv
python -m vmbdiet.cli maaslin \
    --features example/maaslin_demo/data/features.tsv \
    --meta     example/maaslin_demo/data/metadata.tsv \
    --predictors fiber,animal_protein \
    --covariates age \
    --group participant \
    --body-site body_site \
    --out example/maaslin_demo/results
```

## What it does

Two models per feature, each with a **random intercept for participant** — faithful to MaAsLin 3:

- **abundance** — `log(relative abundance) ~ diet + covariates + body_site + (1 | participant)`, fit on
  the samples where the feature is present (statsmodels mixed LM; OLS fallback if the random-effect
  variance is unidentifiable).
- **prevalence** — `present(0/1) ~ same fixed effects`, as a logistic **GEE** with an exchangeable
  within-participant correlation (clustering-aware, population-averaged).

## Planted biology (what the demo should recover)

| feature | planted effect | recovered (`maaslin_py_results.tsv`) |
|---|---|---|
| `PWY_glycogen_degradation` | abundance ↑ with **fiber** | abundance: effect **+0.33**, q ≈ 2e-5 ✅ |
| `SCFA_butyrate_synth` | abundance ↑ with **animal_protein** | abundance: effect **+0.42**, q ≈ 3e-7 ✅ |
| `GUS_beta_glucuronidase` | **prevalence** ↑ with **fiber** | prevalence: effect **+1.67**, q ≈ 2e-12 ✅ |
| `PWY_participant_null` | none (participant-driven) | not significant for diet ✅ |
| `PWY_flat` / `PWY_noise` | none | nominal p<0.05 but **do not survive BH-FDR** (q>0.05) ✅ |

The feature names are mechanistic stand-ins for the Aim-3 estrobolome/SCFA targets (β-glucuronidase,
glycogen/starch degradation, SCFA biosynthesis), so the exact same call runs on real HUMAnN 3 pathway or
eggNOG gene-family tables — just swap the feature table.

## The companion R script

`results/run_maaslin3.R` is the **definitive** fit for the cluster (R + Bioconductor `maaslin3`),
generated from the same arguments as the Python run so the two are directly comparable. The Python
engine is for developing, validating, and demoing the pipeline on a laptop *before* the real tables (or
controlled-access cohorts) are in hand; the R script is what you run once they are.

> Why the random effect matters: `test_random_effect_is_more_conservative_than_naive` in
> `tests/test_maaslin.py` shows that treating the repeated samples as independent (`--group` unset)
> understates p-values for a participant-driven feature — the classic pseudo-replication false positive
> this design has to control.
