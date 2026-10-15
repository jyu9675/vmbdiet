"""Synthetic repeated-measures tests for the MaAsLin 3 bridge.

Builds a gut-vagina-axis-shaped dataset — each participant gives rectal + vaginal samples at two
timepoints — with a planted diet->feature abundance signal and a planted diet->prevalence signal, then
confirms (1) the random-effects engine recovers both, and (2) ignoring the random intercept inflates
significance (the reason MaAsLin 3's mixed model is needed here)."""
import numpy as np
import pandas as pd
import pytest
from vmbdiet import maaslin_py, build_formula, write_maaslin3_script

FEATURES = ["feat_diet_abund", "feat_diet_prev", "feat_clustered_null",
            "feat_flat_hi", "feat_flat_lo", "feat_noise"]


def _repeated_measures(n_part=45, timepoints=2, seed=7):
    """n_part participants x {rectal, vaginal} x timepoints. 'fiber' is a per-participant diet value
    (between-subject, as in an FFQ). Returns (feature_counts, metadata)."""
    rng = np.random.default_rng(seed)
    fiber_by_part = rng.normal(25, 8, n_part)                 # g/day, one per participant
    base_by_part = rng.normal(0, 1.2, n_part)                 # participant random intercept (null feature)
    rows, meta_rows, index = [], [], []
    for p in range(n_part):
        for site in ("rectal", "vaginal"):
            for t in range(timepoints):
                fib = fiber_by_part[p] + rng.normal(0, 1.0)   # tiny within-person visit noise
                # feat_diet_abund: log-abundance rises with fiber (TRUE abundance effect)
                a_signal = np.exp(2.0 + 0.08 * (fib - 25) + rng.normal(0, 0.3))
                # feat_diet_prev: present with prob rising in fiber (TRUE prevalence effect)
                p_present = 1 / (1 + np.exp(-(0.25 * (fib - 25))))
                a_prev = (rng.random() < p_present) * np.exp(1.5 + rng.normal(0, 0.3))
                # feat_clustered_null: NO fiber effect, strong per-participant baseline
                a_null = np.exp(2.0 + 1.5 * base_by_part[p] + rng.normal(0, 0.25))
                a_hi = np.exp(3.0 + rng.normal(0, 0.2))       # always-present filler
                a_lo = (rng.random() < 0.08) * np.exp(1.0)    # rare filler (below min_prevalence)
                a_noise = np.exp(1.5 + rng.normal(0, 0.8))
                rows.append([a_signal, a_prev, a_null, a_hi, a_lo, a_noise])
                meta_rows.append(dict(participant=f"P{p:03d}", body_site=site, timepoint=t,
                                      fiber=fib, age=float(rng.normal(30, 5))))
                index.append(f"P{p:03d}_{site[:3]}_t{t}")
    feat = pd.DataFrame(rows, columns=FEATURES, index=index)
    meta = pd.DataFrame(meta_rows, index=index)
    return feat, meta


def test_recovers_abundance_and_prevalence_signals():
    feat, meta = _repeated_measures()
    res = maaslin_py(feat, meta, predictors=["fiber"], covariates=["age"],
                     group="participant", body_site="body_site")
    assert not res.empty and {"abundance", "prevalence"} <= set(res["model"])

    ab = res[(res.model == "abundance") & (res.feature == "feat_diet_abund") & (res.predictor == "fiber")].iloc[0]
    assert ab["effect"] > 0 and ab["p_value"] < 0.05          # fiber raises abundance, detected

    pv = res[(res.model == "prevalence") & (res.feature == "feat_diet_prev") & (res.predictor == "fiber")].iloc[0]
    assert pv["effect"] > 0 and pv["p_value"] < 0.05          # fiber raises prevalence, detected

    # the always-present filler carries no prevalence variation -> excluded from prevalence model
    assert not ((res.model == "prevalence") & (res.feature == "feat_flat_hi")).any()
    # the rare filler sits below default min_prevalence -> excluded entirely
    assert (res.feature == "feat_flat_lo").sum() == 0
    # multiple-testing column is present per model
    assert res["q_value"].notna().any()


def test_random_effect_is_more_conservative_than_naive():
    """For a feature whose abundance is driven by participant identity (not diet), treating the repeated
    samples as independent (group=None) understates the p-value vs the random-intercept model."""
    feat, meta = _repeated_measures()
    mixed = maaslin_py(feat, meta, predictors=["fiber"], covariates=["age"],
                       group="participant", body_site="body_site", models=("abundance",))
    naive = maaslin_py(feat, meta, predictors=["fiber"], covariates=["age"],
                       group=None, body_site="body_site", models=("abundance",))
    m = mixed[(mixed.feature == "feat_clustered_null") & (mixed.predictor == "fiber")].iloc[0]
    nv = naive[(naive.feature == "feat_clustered_null") & (naive.predictor == "fiber")].iloc[0]
    # naive pseudo-replication -> smaller (anticonservative) p than the clustered-aware fit
    assert nv["p_value"] <= m["p_value"] + 1e-9


def test_formula_and_r_script(tmp_path):
    f = build_formula(["animal_protein", "fiber"], covariates=["age", "race"],
                      body_site="body_site", group="participant")
    assert f == "~ animal_protein + fiber + age + race + body_site + (1 | participant)"
    out = write_maaslin3_script(str(tmp_path / "run_maaslin3.R"),
                                predictors=["animal_protein", "fiber"], covariates=["age", "race"],
                                body_site="body_site", group="participant")
    txt = open(out).read()
    assert "library(maaslin3)" in txt and "(1 | participant)" in txt and 'normalization     = "TSS"' in txt
