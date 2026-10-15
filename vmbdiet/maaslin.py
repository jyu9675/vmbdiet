"""
vmbdiet.maaslin — a MaAsLin 3 bridge for the gut-vagina axis.

The diet-BV literature (and vmbdiet's own `associate` module) tests one sample at a time with
fixed-effects logistic/Spearman models. That is fine for a single vaginal swab per person, but the
FRESH / gut-vagina design is repeated-measures: each participant contributes **two body sites**
(rectal + vaginal) at **up to three timepoints**, so samples from the same host are not independent.
The field-standard tool for this is MaAsLin 3, which (1) fits a per-feature **abundance** model on the
samples where the feature is present and a separate **prevalence** (presence/absence) model, and
(2) absorbs the within-host correlation with a random intercept `(1 | participant)`.

MaAsLin 3 is R-only, so this module does two things:

  * `maaslin_py(...)`  — a faithful *Python* re-implementation of the two-model, random-effects idea
                         (abundance = mixed LM on log relative abundance of present samples; prevalence
                         = logistic GEE clustered by participant). It runs with nothing but statsmodels,
                         so the whole pipeline can be developed and validated on a laptop *before* the
                         real HUMAnN 3 / species tables (or controlled FRESH access) are in hand.
  * `write_maaslin3_script(...)` / `prepare_inputs(...)` — emit MaAsLin 3-format inputs and a
                         ready-to-run R script with the correct formula, for the lab cluster where R +
                         the maaslin3 Bioconductor package live. The Python run and the R run are set up
                         from the *same* arguments, so results are directly comparable.

This is deliberately model-agnostic about the feature table: species relative abundances (Aims 1-2) or
HUMAnN 3 pathway / eggNOG gene-family abundances (Aim 3) go through the same call.
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd

try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    _HAS_SM = True
except Exception:  # pragma: no cover
    _HAS_SM = False


# --------------------------------------------------------------------------- helpers

def _is_categorical(s: pd.Series) -> bool:
    return s.dtype == object or str(s.dtype).startswith("category") or (s.nunique() <= 2 and s.dtype != float)


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR (same implementation as vmbdiet.associate)."""
    p = np.asarray(p, float); n = len(p); order = np.argsort(p)
    q = np.empty(n); prev = 1.0
    for rank, i in enumerate(order[::-1]):
        k = n - rank
        prev = min(prev, p[i] * n / k)
        q[i] = prev
    return q


def tss(feature_table: pd.DataFrame) -> pd.DataFrame:
    """Total-sum-scaling: each sample (row) sums to 1 across features. MaAsLin 3's default."""
    s = feature_table.sum(axis=1)
    return feature_table.div(s.where(s > 0, np.nan), axis=0).fillna(0.0)


def _safe_design(metadata: pd.DataFrame, predictors, covariates, body_site, group, standardize):
    """Rename every model column to a safe token (p0, c0, bs, grp) so the patsy formula never trips
    over spaces/dots in names like 'Energy kcal' or 'animal_protein_density'. Returns (df, formula_rhs,
    token->label map, predictor token list)."""
    cols, rename, labels = {}, {}, {}
    pred_tokens = []
    for i, c in enumerate(predictors):
        tok = f"p{i}"; rename[c] = tok; labels[tok] = c; pred_tokens.append(tok)
    for i, c in enumerate(covariates):
        tok = f"c{i}"; rename[c] = tok; labels[tok] = c
    keep = list(predictors) + list(covariates)
    if body_site:
        rename[body_site] = "bs"; labels["bs"] = body_site; keep.append(body_site)
    if group:
        rename[group] = "grp"; keep.append(group)
    d = metadata[keep].rename(columns=rename).copy()

    terms = []
    for tok in pred_tokens + ([f"c{i}" for i in range(len(covariates))]):
        col = d[tok]
        if _is_categorical(col):
            terms.append(f"C({tok})")
        else:
            if standardize:
                sd = col.std(ddof=0)
                if sd and sd > 0:
                    d[tok] = (col - col.mean()) / sd
            terms.append(tok)
    if body_site:
        terms.append("C(bs)")
    return d, " + ".join(terms), labels, pred_tokens


# --------------------------------------------------------------------------- Python engine

def maaslin_py(feature_table: pd.DataFrame, metadata: pd.DataFrame, predictors,
               covariates=(), group: str | None = None, body_site: str | None = None,
               normalize: str = "tss", min_prevalence: float = 0.1,
               standardize: bool = True, models=("abundance", "prevalence"),
               abundance_pseudocount: float | None = None) -> pd.DataFrame:
    """A MaAsLin 3-style scan in pure Python.

    For every feature passing `min_prevalence`, fit up to two models and report each PREDICTOR's effect:

      * abundance  — on the samples where the feature is PRESENT: log(relative abundance) ~ predictors +
                     covariates [+ body_site], with a random intercept for `group` (mixed LM). Effect =
                     change in log-abundance per predictor unit (per-SD if standardize).
      * prevalence — on ALL samples: 1(feature present) ~ same fixed effects, logistic GEE with an
                     exchangeable within-`group` correlation (population-averaged, clustering-aware).
                     Effect = log-odds of presence per predictor unit.

    If `group` is None both models drop the random/clustered part (plain OLS / logit) — used to show
    what repeated-measures non-independence does to the p-values.

    Returns a tidy long DataFrame: feature, model, predictor, effect, stderr, p_value, q_value (BH within
    model), n, n_present. Per-feature failures are captured in a `note` column, never raised.
    """
    if not _HAS_SM:
        raise RuntimeError("statsmodels required: pip install statsmodels")

    ft = feature_table.copy()
    idx = ft.index.intersection(metadata.index)
    ft = ft.loc[idx]
    meta = metadata.loc[idx]
    rel = tss(ft) if normalize == "tss" else ft
    n = len(idx)
    pc = abundance_pseudocount if abundance_pseudocount is not None else \
        (rel.values[rel.values > 0].min() / 2 if (rel.values > 0).any() else 1e-6)

    d0, rhs, labels, pred_tokens = _safe_design(meta, list(predictors), list(covariates),
                                                body_site, group, standardize)
    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for feat in rel.columns:
            y = rel[feat]
            present = y > 0
            prev = float(present.mean())
            if prev < min_prevalence:
                continue
            # ---- abundance: mixed LM on present samples ----
            if "abundance" in models:
                rows.append(_fit_one("abundance", feat, np.log(y[present] + pc), d0[present.values],
                                     rhs, labels, pred_tokens, group, int(present.sum()), n, prev))
            # ---- prevalence: logistic GEE on all samples ----
            if "prevalence" in models:
                # a feature present (or absent) in everyone carries no prevalence signal
                if 0 < present.sum() < n:
                    rows.append(_fit_one("prevalence", feat, present.astype(float), d0, rhs, labels,
                                         pred_tokens, group, n, n, prev))
    res = pd.DataFrame([r for sub in rows for r in sub])
    if res.empty:
        return res
    for m in res["model"].unique():
        mask = (res["model"] == m) & res["p_value"].notna()
        if mask.any():
            res.loc[mask, "q_value"] = _bh(res.loc[mask, "p_value"].to_numpy())
    return res.sort_values(["model", "p_value"], na_position="last").reset_index(drop=True)


def _fit_one(model, feat, y, d, rhs, labels, pred_tokens, group, n, n_total, prev):
    """Fit one feature's abundance (mixed LM) or prevalence (logistic GEE) model; return a list of
    one result dict per predictor. All failures are caught and reported as a note row."""
    d = d.copy()
    d["y"] = y.values if hasattr(y, "values") else y
    formula = "y ~ " + rhs if rhs else "y ~ 1"
    base = dict(feature=feat, model=model, n=n, n_present=int(round(prev * n_total)),
                prevalence=round(prev, 3))
    fit_note = ""
    try:
        if model == "abundance":
            if group and d["grp"].nunique() > 1 and len(d) > d["grp"].nunique() + 2:
                try:
                    fit = smf.mixedlm(formula, d, groups=d["grp"]).fit(reml=False, method="lbfgs")
                except Exception:  # random-intercept fit failed to converge (e.g. RE variance ~0)
                    fit = smf.ols(formula, d).fit(); fit_note = "RE failed; OLS fallback"
            else:
                fit = smf.ols(formula, d).fit()
        else:  # prevalence
            if y.sum() in (0, len(y)):
                raise ValueError("no prevalence variation")
            if group and d["grp"].nunique() > 1:
                try:
                    fit = smf.gee(formula, groups=d["grp"], data=d,
                                  family=sm.families.Binomial(),
                                  cov_struct=sm.cov_struct.Exchangeable()).fit()
                except Exception:
                    fit = smf.logit(formula, d).fit(disp=0); fit_note = "GEE failed; logit fallback"
            else:
                fit = smf.logit(formula, d).fit(disp=0)
        out = []
        for tok in pred_tokens:
            # a categorical predictor expands to C(tok)[T.level]; report every level found
            hits = [p for p in fit.params.index if p == tok or p.startswith(f"C({tok})")]
            if not hits:
                out.append({**base, "predictor": labels[tok], "effect": np.nan,
                            "stderr": np.nan, "p_value": np.nan, "note": "term dropped"})
            for name in hits:
                lvl = labels[tok] + name[name.find("[T."):].rstrip("]").replace("[T.", "=") if "[T." in name else labels[tok]
                out.append({**base, "predictor": lvl, "effect": float(fit.params[name]),
                            "stderr": float(fit.bse[name]), "p_value": float(fit.pvalues[name]),
                            "note": fit_note})
        return out
    except Exception as e:  # pragma: no cover - robustness path
        return [{**base, "predictor": labels[t], "effect": np.nan, "stderr": np.nan,
                 "p_value": np.nan, "note": str(e)[:70]} for t in pred_tokens]


# --------------------------------------------------------------------------- R / MaAsLin 3 export

def build_formula(predictors, covariates=(), body_site: str | None = None,
                  group: str | None = None) -> str:
    """The lme4-style formula MaAsLin 3 (and lme4/glmmTMB) expect, e.g.
    `~ animal_protein_density + fiber + age + C(race) + body_site + (1 | participant)`."""
    terms = list(predictors) + list(covariates)
    if body_site:
        terms.append(body_site)
    rhs = " + ".join(terms) if terms else "1"
    if group:
        rhs += f" + (1 | {group})"
    return "~ " + rhs


def prepare_inputs(feature_table: pd.DataFrame, metadata: pd.DataFrame, out_dir: str,
                   normalize_abundance: bool = False):
    """Write MaAsLin 3-format inputs: features.tsv (samples x features) and metadata.tsv (samples x
    variables), aligned on the shared sample index. MaAsLin 3 TSS-normalizes internally, so pass raw
    relative abundances / counts (normalize_abundance stays False unless you pre-normalized)."""
    import os
    os.makedirs(out_dir, exist_ok=True)
    idx = feature_table.index.intersection(metadata.index)
    ft = tss(feature_table.loc[idx]) if normalize_abundance else feature_table.loc[idx]
    fpath = os.path.join(out_dir, "features.tsv")
    mpath = os.path.join(out_dir, "metadata.tsv")
    ft.to_csv(fpath, sep="\t")
    metadata.loc[idx].to_csv(mpath, sep="\t")
    return fpath, mpath


def write_maaslin3_script(out_path: str, predictors, covariates=(), body_site: str | None = None,
                          group: str | None = None, features_tsv: str = "features.tsv",
                          metadata_tsv: str = "metadata.tsv", output_dir: str = "maaslin3_out",
                          normalization: str = "TSS", transform: str = "LOG",
                          correction: str = "BH") -> str:
    """Emit a ready-to-run MaAsLin 3 R script wired with the repeated-measures formula. Run it on the
    cluster (R + Bioconductor `maaslin3`) once the real feature/metadata TSVs are in place."""
    formula = build_formula(predictors, covariates, body_site, group)
    script = f'''#!/usr/bin/env Rscript
# Auto-generated by vmbdiet.maaslin - MaAsLin 3 model for the gut-vagina axis.
# install once: BiocManager::install("maaslin3")
suppressMessages(library(maaslin3))

features <- read.delim("{features_tsv}", row.names = 1, check.names = FALSE)
metadata <- read.delim("{metadata_tsv}", row.names = 1, check.names = FALSE)

# Repeated measures: each participant contributes rectal + vaginal samples across timepoints,
# so the random intercept (1 | {group}) absorbs the within-host correlation.
fit <- maaslin3(
    input_data        = features,
    input_metadata    = metadata,
    output            = "{output_dir}",
    formula           = {formula},
    normalization     = "{normalization}",
    transform         = "{transform}",
    correction        = "{correction}",
    augment           = TRUE,
    standardize       = TRUE,                        # per-SD effects, comparable across diet variables
    median_comparison_abundance = TRUE,              # MaAsLin 3 abundance null at the median
    cores             = 4
)
# Outputs: {output_dir}/significant_results.tsv  (abundance + prevalence associations, BH-corrected)
'''
    with open(out_path, "w", newline="\n", encoding="utf-8") as fh:
        fh.write(script)
    return out_path
