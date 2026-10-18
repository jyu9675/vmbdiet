"""
vmbdiet.associate — diet ↔ vaginal-microbiome association tests.

Implements the statistical core shared across the diet–BV literature:
  * `logistic_assoc`    — outcome (BV / CST-IV / Lactobacillus-depleted, 0/1) ~ nutrient + covariates,
                          reported as adjusted odds ratio + 95% CI + p (the primary model in
                          Neggers 2007, Tuddenham 2019, and the HEI/CST studies).
  * `nutrient_taxon_corr` — Spearman/Pearson correlation of a nutrient with each taxon's abundance,
                          with Benjamini-Hochberg FDR (the carbohydrate ↔ L. crispatus style result).
  * `mannwhitney_by_group` — nutrient intake compared between Lactobacillus-dominant vs -depleted.
  * `permanova`          — community composition (a distance matrix) ~ a grouping, permutation F-test
                          (the "does diet structure beta-diversity" test).
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import stats

try:
    import statsmodels.formula.api as smf
    _HAS_SM = True
except Exception:  # pragma: no cover
    _HAS_SM = False


def _is_categorical(s: pd.Series) -> bool:
    return s.dtype == object or str(s.dtype).startswith("category") or s.nunique() <= 2 and s.dtype != float


def logistic_assoc(df: pd.DataFrame, outcome: str, predictor: str,
                   covariates=(), standardize: bool = True) -> dict:
    """Adjusted logistic regression: outcome(0/1) ~ predictor + covariates.

    Categorical covariates are wrapped in C(); a continuous predictor is z-scored by default so the
    OR is per-SD (comparable across nutrients). Returns the PREDICTOR's odds ratio, 95% CI, p, n.
    """
    if not _HAS_SM:
        raise RuntimeError("statsmodels required: pip install statsmodels")
    cols = [outcome, predictor] + list(covariates)
    d = df[cols].dropna().copy()
    d[outcome] = d[outcome].astype(float)
    pterm = predictor
    if standardize and np.issubdtype(d[predictor].dtype, np.number):
        sd = d[predictor].std(ddof=0)
        if sd > 0:
            d["_z"] = (d[predictor] - d[predictor].mean()) / sd
            pterm = "_z"
    terms = [pterm]
    for c in covariates:
        # wrap in Q() so covariate names with spaces/dots (e.g. "Energy kcal") are valid in the formula
        terms.append(f"C(Q('{c}'))" if _is_categorical(d[c]) else f"Q('{c}')")
    formula = f"Q('{outcome}') ~ " + " + ".join(terms)
    res = smf.logit(formula, data=d).fit(disp=0)
    name = pterm if pterm in res.params.index else predictor
    coef = res.params[name]; ci = res.conf_int().loc[name]
    return dict(outcome=outcome, predictor=predictor, n=int(len(d)),
                odds_ratio=float(np.exp(coef)), ci_low=float(np.exp(ci[0])), ci_high=float(np.exp(ci[1])),
                p_value=float(res.pvalues[name]), per="SD" if pterm == "_z" else "unit",
                covariates=list(covariates))


def logistic_scan(df, outcome, predictors, covariates=()) -> pd.DataFrame:
    """Run `logistic_assoc` for many nutrients; BH-FDR across them. Sorted by p."""
    out = []
    for p in predictors:
        try:
            out.append(logistic_assoc(df, outcome, p, covariates))
        except Exception as e:  # keep going; report the failure
            out.append(dict(outcome=outcome, predictor=p, n=0, odds_ratio=np.nan,
                            ci_low=np.nan, ci_high=np.nan, p_value=np.nan, note=str(e)[:60]))
    r = pd.DataFrame(out)
    ok = r["p_value"].notna()
    r.loc[ok, "q_value"] = _bh(r.loc[ok, "p_value"].to_numpy())
    return r.sort_values("p_value", na_position="last").reset_index(drop=True)


def nutrient_taxon_corr(diet: pd.Series, abund: pd.DataFrame, method: str = "spearman") -> pd.DataFrame:
    """Correlate one nutrient (diet, indexed by sample) with each taxon's abundance. BH-FDR."""
    rel = abund.div(abund.sum(axis=1), axis=0).fillna(0.0)
    idx = diet.dropna().index.intersection(rel.index)
    x = diet.loc[idx]
    fn = stats.spearmanr if method == "spearman" else stats.pearsonr
    rows = []
    for t in rel.columns:
        y = rel.loc[idx, t]
        if y.nunique() < 3:
            continue
        r, p = fn(x, y)
        rows.append(dict(taxon=t, r=float(r), p_value=float(p), n=len(idx)))
    df = pd.DataFrame(rows)
    if not df.empty:
        df["q_value"] = _bh(df["p_value"].to_numpy())
        df = df.sort_values("p_value").reset_index(drop=True)
    return df


def mannwhitney_by_group(values: pd.Series, groups: pd.Series) -> dict:
    """Mann-Whitney U of a nutrient between two groups (e.g. Lactobacillus dominant vs depleted)."""
    idx = values.dropna().index.intersection(groups.dropna().index)
    v = values.loc[idx]; g = groups.loc[idx]
    levels = sorted(g.unique())
    if len(levels) != 2:
        raise ValueError(f"need exactly 2 groups, got {levels}")
    a = v[g == levels[0]]; b = v[g == levels[1]]
    U, p = stats.mannwhitneyu(a, b, alternative="two-sided")
    return dict(group_a=str(levels[0]), group_b=str(levels[1]), n_a=len(a), n_b=len(b),
                median_a=float(a.median()), median_b=float(b.median()), U=float(U), p_value=float(p))


def permanova(dist: pd.DataFrame, grouping: pd.Series, permutations: int = 999, seed: int = 0) -> dict:
    """PERMANOVA (Anderson 2001): is community composition different across groups?

    dist: square symmetric distance DataFrame (e.g. Bray-Curtis); grouping: sample -> label.
    Returns pseudo-F and a permutation p-value.
    """
    idx = [s for s in dist.index if s in grouping.dropna().index]
    D = dist.loc[idx, idx].to_numpy(); g = grouping.loc[idx].to_numpy()
    n = len(idx); labels = np.unique(g)
    if n < 3 or len(labels) < 2:
        raise ValueError("need >=3 samples and >=2 groups")
    D2 = D ** 2
    SST = D2[np.triu_indices(n, 1)].sum() / n

    def ssw(gl):
        s = 0.0
        for lab in labels:
            m = gl == lab; k = int(m.sum())
            if k > 1:
                s += D2[np.ix_(m, m)][np.triu_indices(k, 1)].sum() / k
        return s
    a = len(labels)
    SSW = ssw(g); SSA = SST - SSW
    F = (SSA / (a - 1)) / (SSW / (n - a)) if SSW > 0 else np.nan
    rng = np.random.default_rng(seed); ge = 1
    for _ in range(permutations):
        gp = rng.permutation(g); sw = ssw(gp)
        Fp = ((SST - sw) / (a - 1)) / (sw / (n - a)) if sw > 0 else np.nan
        if np.isfinite(Fp) and Fp >= F:
            ge += 1
    return dict(pseudo_F=float(F), p_value=ge / (permutations + 1), n=n,
                groups={str(l): int((g == l).sum()) for l in labels}, permutations=permutations)


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR."""
    p = np.asarray(p, float); n = len(p); order = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for rank, i in enumerate(order[::-1]):
        k = n - rank
        prev = min(prev, p[i] * n / k)
        q[i] = prev
    return q


def or_to_logor(odds_ratio: float, ci_low: float, ci_high: float):
    """Convert a reported OR + 95% CI to (log-OR, SE) for meta-analysis."""
    return float(np.log(odds_ratio)), float((np.log(ci_high) - np.log(ci_low)) / (2 * 1.959964))


def meta_analyze(estimates, method: str = "fixed") -> dict:
    """Inverse-variance meta-analysis of per-study effects on the log scale (e.g. log-ORs).

    estimates: list of (beta, se). method: 'fixed' or 'DL' (DerSimonian-Laird random effects).
    Returns the pooled beta + SE + 95% CI (log scale), z, p, and heterogeneity Q / I^2 / tau^2 / k.
    Exponentiate `beta` and `ci` for an odds ratio. The field's standard tool for combining the same
    association measured in several cohorts (here: diet -> vaginal CST across public cohorts).
    """
    b = np.array([e[0] for e in estimates], float)
    se = np.array([e[1] for e in estimates], float)
    if len(b) < 2:
        raise ValueError("need >=2 studies to meta-analyze")
    w = 1.0 / se ** 2
    b_fixed = float((w * b).sum() / w.sum())
    k = len(b); df = k - 1
    Q = float((w * (b - b_fixed) ** 2).sum())
    I2 = max(0.0, (Q - df) / Q) * 100 if Q > 0 else 0.0
    denom = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (Q - df) / denom) if (df > 0 and denom > 0) else 0.0
    if method == "DL":
        wr = 1.0 / (se ** 2 + tau2)
        beta = float((wr * b).sum() / wr.sum()); se_p = float((1.0 / wr.sum()) ** 0.5)
    elif method == "fixed":
        beta = b_fixed; se_p = float((1.0 / w.sum()) ** 0.5)
    else:
        raise ValueError("method must be 'fixed' or 'DL'")
    z = beta / se_p
    p = float(2 * stats.norm.sf(abs(z)))
    ci = (beta - 1.959964 * se_p, beta + 1.959964 * se_p)
    return dict(beta=beta, se=se_p, ci=(float(ci[0]), float(ci[1])), z=float(z), p_value=p,
                Q=Q, I2=float(I2), tau2=float(tau2), k=k, method=method)
