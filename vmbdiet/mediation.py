"""
vmbdiet.mediation — does diet shape the vaginal community *through* the gut?

MaAsLin (vmbdiet.maaslin) tests the marginal links — diet<->gut function, diet<->vaginal CST — but the
gut-vagina axis is a *mediation* hypothesis:

        diet (X)  ->  gut microbial function (M)  ->  vaginal community (Y)
                         e.g. gut beta-glucuronidase / glycogen / SCFA capacity   e.g. CST-IV, L. crispatus

Mediation decomposes the total effect of diet on the vaginal outcome into the part that flows **through**
the gut mediator (natural indirect effect, NIE) and the part that does not (natural direct effect, NDE),
and reports the **proportion mediated**. That is the quantity that actually answers "whether diet shapes
the vaginal microbiome through the gut", which no marginal association can.

Implementation (counterfactual / potential-outcomes framework, Imai-Keele-Tingley 2010; VanderWeele):
  * mediator model  M ~ X + covariates            (linear, OLS)
  * outcome model   Y ~ X + M + covariates        (linear if Y continuous; logistic if Y binary)
  * continuous Y    -> closed-form NIE = a*b, NDE = c'  (exact for linear-linear).
  * binary Y        -> simulation estimator on the **risk-difference scale**: average over the sample of
                       P(Y | X=x1, M(x1)) - P(Y | X=x1, M(x0))  (NIE) and
                       P(Y | X=x1, M(x0)) - P(Y | X=x0, M(x0))  (NDE),
                       with the mediator drawn under each exposure using common random numbers. This is
                       correct when the outcome model is logistic (the naive product-of-coefficients is
                       not, because of non-collapsibility).
Exposure contrast defaults to **+1 SD** of X (X is z-scored), so effects read "per 1-SD increase in diet".
Confidence intervals come from a **cluster bootstrap that resamples participants**, which respects the
repeated-measures non-independence (two body sites x multiple timepoints per host) without a full
multilevel-mediation model.

CAUSAL ASSUMPTIONS — read before interpreting (these are strong and mostly untestable here):
  1. **Temporal order X -> M -> Y.** Cross-sectional data cannot establish it; a mediation estimate on
     same-visit samples is associational. Use diet/function at an earlier visit and CST at a later one
     (pass the right rows via `df`) for a defensible longitudinal version.
  2. No unmeasured confounding of X-Y, X-M, or M-Y, and no M-Y confounder affected by X. Adjust for every
     plausible confounder via `covariates`; residual confounding biases NIE/NDE.
  3. The mediator is measured without substantial error and at the gut site; the outcome at the vaginal
     site. Cross-site is fine; mis-specification of the mediator model is not.
This module computes the estimand; it does not make the assumptions true. Report it as "consistent with
mediation under these assumptions", never as established causation.
"""
from __future__ import annotations
import numpy as np
import pandas as pd

try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from scipy.special import expit
    _HAS_SM = True
except Exception:  # pragma: no cover
    _HAS_SM = False


def _is_categorical(s: pd.Series) -> bool:
    return s.dtype == object or str(s.dtype).startswith("category") or (s.nunique() <= 2 and s.dtype != float)


def _design(df, exposure, mediator, outcome, covariates, group):
    """Rename to safe tokens; z-score the exposure (per-SD contrast). Returns (d, cov_terms, grp_col)."""
    ren = {exposure: "x", mediator: "m", outcome: "y"}
    keep = [exposure, mediator, outcome]
    cov_terms = []
    for i, c in enumerate(covariates):
        ren[c] = f"c{i}"; keep.append(c)
    if group:
        ren[group] = "grp"; keep.append(group)
    d = df[keep].dropna().rename(columns=ren).copy()
    sd = d["x"].std(ddof=0)
    if sd and sd > 0:
        d["x"] = (d["x"] - d["x"].mean()) / sd
    d["y"] = d["y"].astype(float)
    for i, c in enumerate(covariates):
        tok = f"c{i}"
        cov_terms.append(f"C({tok})" if _is_categorical(d[tok]) else tok)
    return d, cov_terms, ("grp" if group else None)


def _point(d, cov_terms, binary, x0, x1, sims, rng):
    """One mediation point estimate on `d`. Returns (total, nde, nie, a, b, cprime)."""
    rhs_cov = (" + " + " + ".join(cov_terms)) if cov_terms else ""
    med = smf.ols(f"m ~ x{rhs_cov}", d).fit()
    a = float(med.params["x"])
    if binary:
        out = smf.logit(f"y ~ x + m{rhs_cov}", d).fit(disp=0)
        b = float(out.params["m"]); cprime = float(out.params["x"])
        # precompute the covariate linear predictor ONCE; the sims loop is then pure numpy
        # (statsmodels .predict re-runs patsy on every call, which is what made this slow).
        names = list(out.model.exog_names); beta = out.params.to_numpy()
        xi, mi = names.index("x"), names.index("m")
        bx, bm = beta[xi], beta[mi]
        full_lp = out.model.exog @ beta
        cov_lp = full_lp - bx * out.model.exog[:, xi] - bm * out.model.exog[:, mi]
        sigma = float(np.sqrt(med.scale))
        mu0 = med.predict(d.assign(x=x0)).to_numpy()
        mu1 = med.predict(d.assign(x=x1)).to_numpy()
        nde_acc = nie_acc = 0.0
        for _ in range(sims):
            eps = rng.normal(0, sigma, len(d))
            m_x0 = mu0 + eps; m_x1 = mu1 + eps               # common random numbers
            p_x1_m0 = expit(cov_lp + bx * x1 + bm * m_x0)
            p_x0_m0 = expit(cov_lp + bx * x0 + bm * m_x0)
            p_x1_m1 = expit(cov_lp + bx * x1 + bm * m_x1)
            nde_acc += float((p_x1_m0 - p_x0_m0).mean())
            nie_acc += float((p_x1_m1 - p_x1_m0).mean())
        nde = nde_acc / sims; nie = nie_acc / sims
        total = nde + nie
    else:
        out = smf.ols(f"y ~ x + m{rhs_cov}", d).fit()
        b = float(out.params["m"]); cprime = float(out.params["x"])
        nie = a * b * (x1 - x0)
        nde = cprime * (x1 - x0)
        total = nie + nde
    return total, nde, nie, a, b, cprime


def mediate(df: pd.DataFrame, exposure: str, mediator: str, outcome: str,
            covariates=(), group: str | None = None, outcome_type: str = "auto",
            contrast: float = 1.0, n_boot: int = 500, sims: int = 100, seed: int = 0) -> dict:
    """Estimate the mediation of `exposure` -> `mediator` -> `outcome`.

    outcome_type: 'auto' (binary if <=2 unique values else continuous), 'binary', or 'continuous'.
    contrast: exposure change in SD units (default +1 SD).
    group: participant column -> cluster bootstrap resamples participants for the CIs.
    Returns effects (total / direct NDE / indirect NIE / proportion_mediated) each with a 95% CI, the
    path coefficients a / b / c', the scale, n, and the assumptions reminder.
    """
    if not _HAS_SM:
        raise RuntimeError("statsmodels required: pip install statsmodels")
    d, cov_terms, grp = _design(df, exposure, mediator, outcome, covariates, group)
    binary = (d["y"].nunique() <= 2) if outcome_type == "auto" else (outcome_type == "binary")
    x0, x1 = 0.0, float(contrast)
    rng = np.random.default_rng(seed)

    total, nde, nie, a, b, cprime = _point(d, cov_terms, binary, x0, x1, sims, rng)
    prop = nie / total if total != 0 else np.nan

    # cluster bootstrap (resample participants if grouped, else rows)
    boots = {"total": [], "nde": [], "nie": [], "prop": []}
    if grp:
        groups = d[grp].unique()
        get = lambda rb: pd.concat([d[d[grp] == g] for g in rng.choice(groups, len(groups), replace=True)])
    else:
        get = lambda rb: d.sample(len(d), replace=True, random_state=rb)
    ok = 0
    for bi in range(n_boot):
        try:
            db = get(bi)
            t, nd, ni, *_ = _point(db, cov_terms, binary, x0, x1, sims, rng)
            p = ni / t if t != 0 else np.nan
            boots["total"].append(t); boots["nde"].append(nd); boots["nie"].append(ni); boots["prop"].append(p)
            ok += 1
        except Exception:
            continue

    def ci(key, est):
        arr = np.array([v for v in boots[key] if np.isfinite(v)])
        if len(arr) < max(20, 0.5 * ok):
            return (float("nan"), float("nan"))
        return (float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5)))

    scale = "risk difference (probability)" if binary else f"units of {outcome}"
    return {
        "exposure": exposure, "mediator": mediator, "outcome": outcome,
        "outcome_type": "binary" if binary else "continuous", "scale": scale,
        "contrast": f"+{contrast} SD of {exposure}", "n": int(len(d)),
        "n_boot_ok": ok, "clustered_by": group,
        "total_effect": total, "total_ci": ci("total", total),
        "direct_effect_nde": nde, "nde_ci": ci("nde", nde),
        "indirect_effect_nie": nie, "nie_ci": ci("nie", nie),
        "proportion_mediated": prop, "prop_ci": ci("prop", prop),
        "path_a_x_to_m": a, "path_b_m_to_y": b, "path_cprime_direct": cprime,
        "assumptions": "Estimand only; valid under temporal order X->M->Y and no unmeasured "
                       "X-Y/X-M/M-Y confounding. Report as 'consistent with mediation', not causation.",
    }
