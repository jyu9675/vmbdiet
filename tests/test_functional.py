"""Tests for the Aim-3 functional layer.

Parser/matcher is validated on a REAL (trimmed) HUMAnN 3 gene-family table from biobakery's demo data;
the dictionary and the functional -> MaAsLin end-to-end are validated on a synthetic EC-level table with
planted diet->function signals.
"""
import re
import numpy as np
import pandas as pd
import pytest
from pathlib import Path
from vmbdiet import (read_humann_table, detect_level, match_features, build_axis_feature_table,
                     FUNCTIONAL_FEATURES, CATEGORIES, maaslin_py, write_cmd_export_R)

REAL_EXCERPT = Path(__file__).parent.parent / "example" / "functional_demo" / "data" / "demo_genefamilies_excerpt.tsv"


def test_dictionary_integrity():
    assert set(CATEGORIES) == {"estrobolome", "glycogen_starch", "scfa"}
    for cat, spec in FUNCTIONAL_FEATURES.items():
        assert spec["label"] and spec["note"]
        for ec in spec["ecs"]:
            assert re.fullmatch(r"\d+\.\d+\.\d+\.\d+", ec), f"bad EC in {cat}: {ec}"
        for ko in spec["kos"]:
            assert re.fullmatch(r"K\d{5}", ko), f"bad KO in {cat}: {ko}"
        for pw in spec["pathways"]:
            assert "PWY" in pw, f"bad pathway in {cat}: {pw}"
    # the verified, grounded identifiers are present
    assert "3.2.1.31" in FUNCTIONAL_FEATURES["estrobolome"]["ecs"]
    assert "K01195" in FUNCTIONAL_FEATURES["estrobolome"]["kos"]
    assert "PWY-5676" in FUNCTIONAL_FEATURES["scfa"]["pathways"]
    assert "GLYCOCAT-PWY" in FUNCTIONAL_FEATURES["glycogen_starch"]["pathways"]


def test_detect_level():
    assert detect_level(["K01195", "K00925", "K00929"]) == "ko"
    assert detect_level(["3.2.1.31", "2.7.2.7"]) == "ec"
    assert detect_level(["PWY-5676", "GLYCOCAT-PWY", "CENTFERM-PWY"]) == "pathway"
    assert detect_level(["UniRef50_A6L108", "UniRef50_R5FJB9"]) == "uniref"


@pytest.mark.skipif(not REAL_EXCERPT.exists(), reason="real HUMAnN excerpt fixture missing")
def test_parse_and_match_real_humann():
    t = read_humann_table(REAL_EXCERPT)                      # default stratified='total'
    # stratification + special rows are dropped from community totals
    assert not t.index.str.contains(r"\|").any()
    assert "UNMAPPED" not in t["feature_id"].values
    assert detect_level(t["feature_id"]) == "uniref"

    m = match_features(t)
    cats = set(m["category"])
    assert "scfa" in cats and "glycogen_starch" in cats     # acetate kinase / amylase are really in the demo
    assert "estrobolome" not in cats                         # no beta-glucuronidase in this Bacteroides demo
    names = " ".join(str(x).lower() for x in m["feature"])
    assert "acetate kinase" in names
    assert ("amylase" in names) or ("glucosidase" in names)
    assert (m["matched_by"] == "name").all()                 # UniRef table -> name matching


def _synthetic_ec_table(n_part=40, seed=3):
    """features x samples EC table; alpha-amylase (glycogen_starch) rises with fiber."""
    rng = np.random.default_rng(seed)
    fiber = rng.normal(25, 8, n_part)
    samples, meta = [], []
    for p in range(n_part):
        for site in ("rectal", "vaginal"):
            samples.append(f"P{p:03d}_{site[:3]}")
            meta.append(dict(participant=f"P{p:03d}", body_site=site,
                             fiber=round(fiber[p] + rng.normal(0, 1), 2)))
    meta = pd.DataFrame(meta, index=samples)
    rows = {
        "3.2.1.1: alpha-amylase": [np.exp(2 + 0.06 * (meta.loc[s, "fiber"] - 25) + rng.normal(0, .3)) for s in samples],
        "3.2.1.31: beta-glucuronidase": [np.exp(2 + rng.normal(0, .3)) for _ in samples],
        "2.7.7.7: DNA polymerase": [np.exp(2 + rng.normal(0, .3)) for _ in samples],
    }
    feat = pd.DataFrame(rows, index=samples).T
    feat.loc["UNMAPPED"] = rng.uniform(1e3, 2e3, len(samples))
    feat.index.name = "# Gene Family"
    return feat, meta


def test_axis_table_and_maaslin_endtoend():
    feat, meta = _synthetic_ec_table()
    t = read_humann_table(feat)
    axis = build_axis_feature_table(t, aggregate="sum")
    assert list(axis.index) == list(feat.columns)            # samples x category
    assert "glycogen_starch" in axis.columns and "estrobolome" in axis.columns

    # feed the axis scores into the random-effects engine (absolute scores -> normalize='none')
    res = maaslin_py(axis, meta, predictors=["fiber"], group="participant",
                     body_site="body_site", normalize="none", models=("abundance",))
    gs = res[(res.feature == "glycogen_starch") & (res.predictor == "fiber")].iloc[0]
    assert gs["effect"] > 0 and gs["p_value"] < 0.05         # planted fiber->glycogen signal recovered
    eb = res[(res.feature == "estrobolome") & (res.predictor == "fiber")].iloc[0]
    assert eb["p_value"] > 0.05                              # no planted estrobolome->fiber effect


def test_cmd_export_r(tmp_path):
    out = write_cmd_export_R(str(tmp_path / "cmd.R"), study="AsnicarF_2021", data_type="pathway_abundance")
    txt = open(out).read()
    assert "curatedMetagenomicData" in txt and "AsnicarF_2021.pathway_abundance" in txt
