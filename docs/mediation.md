# Mediation — testing the axis's causal story

`vmbdiet.mediation` answers the gut–vagina-axis proposal's central question: not *is diet associated
with the vaginal community* (that's [`maaslin`](maaslin-bridge.md)), but *does diet act on the vaginal
community **through** the gut*.

```
diet (X)  ->  gut microbial function (M)  ->  vaginal community (Y)
```

It decomposes the total effect of diet on the vaginal outcome into:

- **NIE** (natural indirect effect) — the part that flows **through** the gut mediator;
- **NDE** (natural direct effect) — the part that does not;
- **proportion mediated** = NIE / total.

## Estimator

Counterfactual / potential-outcomes framework (Imai–Keele–Tingley 2010; VanderWeele):

- mediator model `M ~ X + covariates` (OLS);
- outcome model `Y ~ X + M + covariates` — linear if Y is continuous (`lacto_frac`), **logistic** if Y is
  binary (`cst_IV`);
- **continuous Y**: closed-form NIE = a·b, NDE = c' (exact for linear–linear);
- **binary Y**: a simulation estimator on the **risk-difference scale** (average change in P(CST-IV)),
  because the naive product-of-coefficients is invalid under a logistic outcome (non-collapsibility).
  The mediator is drawn under each exposure with common random numbers; the covariate linear predictor is
  precomputed so the simulation is pure numpy (fast).

The exposure defaults to a **+1 SD** contrast (X is z-scored), so effects read "per 1-SD increase in
diet", comparable to the per-SD MaAsLin effects.

## Repeated measures

Confidence intervals come from a **cluster bootstrap that resamples participants** (`--group
participant`), which respects the two-body-site × multi-timepoint non-independence without a full
multilevel-mediation model. This is the pragmatic, defensible choice; a formal multilevel mediation is a
future extension.

## Causal assumptions — do not skip

The module computes the estimand; it cannot make these true. Interpret results as "consistent with
mediation under these assumptions", never as established causation.

1. **Temporal order X → M → Y.** Cross-sectional (same-visit) data cannot establish it. For a defensible
   longitudinal version, feed diet/function from an earlier visit and CST from a later one.
2. **No unmeasured confounding** of X–Y, X–M, or M–Y, **and no M–Y confounder affected by X.** Adjust for
   every plausible confounder via `covariates` (sexual activity, antibiotics, hormones, BV/STI, BMI,
   age, race); residual confounding biases NIE/NDE.
3. The mediator is measured well (gut site) and the outcome at the vaginal site. Cross-site is fine;
   mediator-model mis-specification is not.

## Why it matters for the proposal

Two cohorts can show diet↔CST *and* diet↔gut-function associations and still leave the mechanism
untested — the associations could be independent. The proportion mediated is the number that directly
supports or undercuts "diet shapes the vaginal ecosystem through the gut", and it is what distinguishes
this project from a pair of marginal scans. See [example/mediation_demo](../example/mediation_demo).
