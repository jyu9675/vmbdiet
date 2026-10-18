# Cross-cohort meta-analysis — diet → vaginal CST-IV

Pools the diet→CST-IV association across the two public diet↔vaginal 16S cohorts analyzed in this repo
([italy_stressdiet](../italy_stressdiet), [melody_ibd](../melody_ibd)) with `vmbdiet.meta_analyze`
(inverse-variance). The point: either cohort alone is limited — MELODY (n=47) is underpowered — but pooling
recovers a consistent signal.

## Run

```bash
python example/cross_cohort_meta/meta.py      # writes meta_results.tsv
```

## Finding (n = 113 + 47 = 160 women, per SD, energy-adjusted)

| diet feature | Italy OR | MELODY OR | **pooled OR [95% CI]** | p | I² |
|---|---|---|---|---|---|
| **animal protein** | 1.66 | 1.30 | **1.54 [1.08, 2.20]** | **0.017** | 0% |
| vegetable protein | 0.68 | 0.81 | 0.72 [0.49, 1.05] | 0.091 | 0% |
| carbohydrate | 1.16 | 0.80 | 1.05 [0.75, 1.48] | 0.78 | 0% |
| fiber | 0.83 | 1.05 | 0.89 [0.63, 1.27] | 0.52 | 0% |
| alcohol | 1.12 | 1.15 | 1.13 [0.82, 1.56] | 0.46 | 0% |
| saturated fat | 0.95 | 0.87 | 0.93 [0.66, 1.31] | 0.66 | 0% |

Higher **animal-protein density → higher odds of the dysbiotic CST-IV state**, significant and consistent
across two independent cohorts (I² = 0%). A suggestive protective trend for vegetable protein (p = 0.09)
gives a coherent animal-vs-plant contrast; fiber, carbohydrate, alcohol and saturated fat show no pooled
signal.

## Caveats (stated up front)

- **Two cohorts only** — a minimal meta-analysis; I² is unstable at k = 2, so "consistent direction" is the
  defensible claim, not a precise heterogeneity estimate.
- Pooled at the **derived CST-IV outcome with standardized (per-SD) exposures**, not raw taxa, because the
  cohorts use **different 16S regions** (pooling features across regions would confound the comparison).
- Observational, energy-adjusted; vaginal 16S only. This covers Aim 1–2 (diet↔community). The mechanistic
  layers — diet → gut function → vaginal community, and rectal↔vaginal strain transfer — need paired
  gut+vaginal+diet+shotgun data, which no public cohort provides.
