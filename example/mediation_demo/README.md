# Mediation demo — does diet act on the vagina *through* the gut?

The capstone of the gut–vagina-axis analysis: Aims 1–2 (`vmbdiet maaslin`) establish the marginal links
diet↔gut-function and diet↔CST, but only **mediation** answers the proposal's actual question —
*whether diet shapes the vaginal microbiome through the gut*.

This demo plants a known partial mediation and recovers it:

```
fiber (X)  ->  gut SCFA capacity (M)  ->  CST-IV (Y)
   fiber raises gut SCFA (a>0);  gut SCFA lowers CST-IV risk (b<0, the mediated path);
   fiber also has a small DIRECT protective effect (c'<0).
```

## Run it

```bash
python example/mediation_demo/make_demo.py        # writes data/axis_table.tsv
python -m vmbdiet.cli mediate \
    --data example/mediation_demo/data/axis_table.tsv \
    --exposure fiber --mediator gut_scfa --outcome cst_IV \
    --covariates age --n-boot 400 \
    --out example/mediation_demo/results/mediation.json
```

## What it finds

| effect | estimate (Δ risk of CST-IV per +1 SD fiber) | 95% CI |
|---|---|---|
| total | −0.26 | [−0.31, −0.21] |
| direct (NDE) | −0.11 | [−0.19, −0.02] |
| **indirect, through gut SCFA (NIE)** | **−0.15** | **[−0.22, −0.09]** |
| proportion mediated | **0.59** | [0.34, 0.91] |

The indirect-effect CI excludes 0 → **~59% of fiber's protection against CST-IV flows through gut SCFA
capacity**, exactly the axis mechanism. On the binary (CST) outcome the effects are on the
**risk-difference** scale via the counterfactual estimator (the naive product-of-coefficients is invalid
for a logistic outcome).

## In a real run

- The **mediator** column is a `vmbdiet functional` axis capacity score computed on the **gut** sample
  (e.g. gut β-glucuronidase, glycogen, or SCFA capacity).
- The **outcome** is a `vmbdiet cst` call on the **vaginal** sample (`cst_IV`, or `lacto_frac` for a
  continuous outcome).
- Pass `--group participant` so the bootstrap resamples participants (the repeated-measures / two-site
  non-independence).

## Read the caveats

Mediation estimates a quantity; it does not make the causal assumptions true. This demo is
cross-sectional, so the temporal order X→M→Y is assumed, not shown. For a defensible longitudinal
version use diet/function from an earlier visit and CST from a later one. Full assumptions in
[docs/mediation.md](../../docs/mediation.md). Report results as "consistent with mediation", never as
established causation.
