"""Synthetic tests: inject a diet->microbiome signal and confirm the tools recover it."""
import numpy as np
import pandas as pd
import pytest
from vmbdiet import (assign_cst, lactobacillus_fraction, logistic_assoc, logistic_scan,
                     nutrient_taxon_corr, mannwhitney_by_group, permanova, bray_curtis, pcoa)

TAXA = ["Lactobacillus_crispatus", "Lactobacillus_iners", "Gardnerella_vaginalis",
        "Prevotella_bivia", "Atopobium_vaginae"]


def _synthetic(n=120, seed=1):
    """Build n samples where higher 'carbohydrate' -> more L. crispatus (CST I),
    and higher 'fat' -> more Gardnerella/CST IV. Returns (abund, meta)."""
    rng = np.random.default_rng(seed)
    carb = rng.normal(250, 60, n)          # g/day
    fat = rng.normal(70, 20, n)
    age = rng.normal(30, 5, n)
    race = rng.choice(["A", "B", "C"], n)
    # latent crispatus propensity rises with carb, falls with fat
    z = 0.02 * (carb - 250) - 0.03 * (fat - 70) + rng.normal(0, 0.5, n)
    rows = []
    for i in range(n):
        if z[i] > 0.3:       # crispatus-dominant
            comp = [rng.uniform(0.7, 0.95), 0.05, 0.02, 0.02, 0.01]
        elif z[i] < -0.3:    # Gardnerella/diverse (CST IV)
            comp = [0.03, 0.05, rng.uniform(0.3, 0.5), rng.uniform(0.2, 0.3), 0.15]
        else:                # iners-dominant
            comp = [0.05, rng.uniform(0.6, 0.85), 0.05, 0.03, 0.02]
        comp = np.array(comp); comp = comp / comp.sum()
        rows.append(comp * rng.integers(5000, 20000))   # as counts
    abund = pd.DataFrame(rows, columns=TAXA, index=[f"S{i:03d}" for i in range(n)])
    meta = pd.DataFrame({"carbohydrate": carb, "fat": fat, "age": age, "race": race},
                        index=abund.index)
    return abund, meta


def test_assign_cst_obvious():
    abund = pd.DataFrame([
        [90, 5, 2, 2, 1],   # crispatus -> I
        [5, 85, 5, 3, 2],   # iners -> III
        [3, 5, 45, 30, 17], # Gardnerella/diverse -> IV
    ], columns=TAXA, index=["a", "b", "c"])
    cst = assign_cst(abund)
    assert cst.loc["a", "cst"] == "I"
    assert cst.loc["b", "cst"] == "III"
    assert cst.loc["c", "cst"] == "IV"
    assert cst.loc["a", "lacto_dominant"] and not cst.loc["c", "lacto_dominant"]


def test_lacto_fraction():
    abund, _ = _synthetic(30)
    lf = lactobacillus_fraction(abund)
    assert (lf >= 0).all() and (lf <= 1.0001).all()


def test_logistic_recovers_fat_bv_signal():
    abund, meta = _synthetic(150)
    cst = assign_cst(abund)
    df = meta.join(cst[["cst", "lacto_dominant"]])
    df["cst_IV"] = (df["cst"] == "IV").astype(int)
    # fat was injected to raise CST-IV risk -> adjusted OR > 1, significant
    res = logistic_assoc(df, "cst_IV", "fat", covariates=["age", "race"])
    assert res["odds_ratio"] > 1.0
    assert res["p_value"] < 0.05
    # carbohydrate protects -> OR < 1
    res2 = logistic_assoc(df, "cst_IV", "carbohydrate", covariates=["age", "race"])
    assert res2["odds_ratio"] < 1.0


def test_scan_fdr_and_corr():
    abund, meta = _synthetic(150)
    cst = assign_cst(abund)
    df = meta.join(cst[["cst"]]); df["cst_IV"] = (df["cst"] == "IV").astype(int)
    scan = logistic_scan(df, "cst_IV", ["fat", "carbohydrate", "age"], ["race"])
    assert "q_value" in scan.columns and len(scan) == 3
    # carbohydrate should correlate positively with L. crispatus abundance
    c = nutrient_taxon_corr(meta["carbohydrate"], abund)
    cris = c[c.taxon == "Lactobacillus_crispatus"].iloc[0]
    assert cris["r"] > 0


def test_permanova_detects_structure():
    abund, meta = _synthetic(120)
    cst = assign_cst(abund)
    D = bray_curtis(abund)
    res = permanova(D, cst["cst"], permutations=199)
    assert res["p_value"] < 0.05 and res["pseudo_F"] > 1
    coords, prop = pcoa(D)
    assert coords.shape[0] == len(abund) and len(prop) >= 1


def test_mannwhitney_group():
    abund, meta = _synthetic(120)
    cst = assign_cst(abund)
    r = mannwhitney_by_group(meta["carbohydrate"], cst["lacto_dominant"].map({True: "dom", False: "dep"}))
    assert set([r["group_a"], r["group_b"]]) == {"dom", "dep"}
