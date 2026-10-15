"""Mediation tests: recover a known linear decomposition, and detect / rule out mediation for a binary
(CST-style) outcome, including the cluster bootstrap for repeated measures."""
import numpy as np
import pandas as pd
import pytest
from vmbdiet import mediate


def test_linear_recovers_known_decomposition():
    rng = np.random.default_rng(0)
    n = 500
    X = rng.normal(0, 1, n)
    M = 0.8 * X + rng.normal(0, 1, n)                 # a = 0.8
    Y = 0.5 * X + 1.2 * M + rng.normal(0, 1, n)       # c' = 0.5, b = 1.2 -> NIE = 0.96, total = 1.46
    df = pd.DataFrame(dict(diet=X, gutfun=M, lacto=Y))
    r = mediate(df, "diet", "gutfun", "lacto", n_boot=120, seed=1)
    assert r["outcome_type"] == "continuous"
    assert abs(r["path_a_x_to_m"] - 0.8) < 0.15 and abs(r["path_b_m_to_y"] - 1.2) < 0.15
    assert abs(r["indirect_effect_nie"] - r["path_a_x_to_m"] * r["path_b_m_to_y"]) < 1e-6  # NIE = a*b
    assert 0.55 < r["proportion_mediated"] < 0.78     # ~0.66
    assert r["nie_ci"][0] > 0                          # indirect effect significant


def test_binary_detects_real_mediation():
    rng = np.random.default_rng(1)
    n = 600
    X = rng.normal(0, 1, n)
    M = 0.9 * X + rng.normal(0, 1, n)
    p = 1 / (1 + np.exp(-(0.3 * X + 1.3 * M)))
    Y = (rng.random(n) < p).astype(int)
    df = pd.DataFrame(dict(diet=X, gutfun=M, cst_IV=Y))
    r = mediate(df, "diet", "gutfun", "cst_IV", n_boot=120, sims=60, seed=2)
    assert r["outcome_type"] == "binary" and "risk difference" in r["scale"]
    assert r["indirect_effect_nie"] > 0
    assert r["nie_ci"][0] > 0                          # CI excludes 0 -> mediation detected
    assert 0.0 < r["proportion_mediated"] <= 1.0


def test_binary_rules_out_absent_mediation():
    rng = np.random.default_rng(2)
    n = 600
    X = rng.normal(0, 1, n)
    M = rng.normal(0, 1, n)                            # mediator independent of diet -> no indirect path
    p = 1 / (1 + np.exp(-(0.3 * X + 1.3 * M)))
    Y = (rng.random(n) < p).astype(int)
    df = pd.DataFrame(dict(diet=X, gutfun=M, cst_IV=Y))
    r = mediate(df, "diet", "gutfun", "cst_IV", n_boot=120, sims=60, seed=3)
    lo, hi = r["nie_ci"]
    assert lo <= 0 <= hi                               # CI covers 0 -> no evidence of mediation


def test_cluster_bootstrap_by_participant():
    rng = np.random.default_rng(3)
    n_part = 60
    rows = []
    for pi in range(n_part):
        base = rng.normal(0, 1)                        # per-participant diet level
        for _ in range(3):                             # 3 repeated samples per participant
            x = base + rng.normal(0, 0.2)
            m = 0.8 * x + rng.normal(0, 1)
            y = 0.4 * x + 1.0 * m + rng.normal(0, 1)
            rows.append(dict(participant=f"P{pi}", diet=x, gutfun=m, lacto=y))
    df = pd.DataFrame(rows)
    r = mediate(df, "diet", "gutfun", "lacto", group="participant", n_boot=120, seed=4)
    assert r["clustered_by"] == "participant" and r["n_boot_ok"] > 50
    assert r["nie_ci"][0] < r["indirect_effect_nie"] < r["nie_ci"][1]  # estimate inside its CI
