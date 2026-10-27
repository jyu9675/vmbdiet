#!/usr/bin/env python3
"""Cross-cohort, reproducibility-first association of gut functional CAPACITY with a cardiometabolic
outcome -- Route A of docs/diet-gut-health-landscape.md.

For each independent cohort it (1) builds the gut-vagina-axis functional capacities (glycogen/starch,
SCFA; estrobolome needs the gene-family export) from a HUMAnN pathway table, (2) runs an adjusted
per-SD logistic association of each capacity with the outcome (e.g. Type 2 Diabetes), then (3) pools
the per-cohort log-odds across cohorts with vmbdiet.meta_analyze (DerSimonian-Laird random effects)
and reports heterogeneity (I^2) and whether the effect REPLICATES (same direction in every cohort AND
a significant pooled estimate). The point is rigor -- an association only "counts" if it reproduces
across independent cohorts -- not a new biological claim.

This is an ASSOCIATION, not a mediation: cMD carries no detailed diet, and its only diet/lifestyle
proxy (`westernized`) never co-varies with a cardiometabolic outcome in the same cohort (every cMD
cardiometabolic cohort is westernized=Yes). So this harness tests function <-> outcome reproducibility;
the diet -> function -> health mediation needs a cohort carrying diet + function + outcome together
(MetaCardis / American Gut / PREDICT -- see the landscape doc).

Get the real inputs with ONE R step per cohort (no further R afterwards), via write_cmd_export_R():
  se <- curatedMetagenomicData("2021-10-14.<Study>.pathway_abundance", dryrun=FALSE, counts=FALSE)[[1]]
  write.table(as.data.frame(assay(se)),   "cmd_<Study>_pathway_abundance.tsv", sep="\t", col.names=NA)
  write.table(as.data.frame(colData(se)), "cmd_<Study>_metadata.tsv",          sep="\t", col.names=NA)
Real T2D cohorts in cMD: KarlssonFH_2013, QinJ_2012, SankaranarayananK_2015 (disease column carries
"Type 2 Diabetes Mellitus" / "Healthy"). Then, e.g.:
  python analyze.py \
    --cohort Karlsson cmd_KarlssonFH_2013_pathway_abundance.tsv cmd_KarlssonFH_2013_metadata.tsv \
    --cohort Qin      cmd_QinJ_2012_pathway_abundance.tsv        cmd_QinJ_2012_metadata.tsv \
    --outcome disease --case "Type 2 Diabetes Mellitus" --control Healthy --covars age,sex,BMI

With no --cohort arguments it generates and runs the committed synthetic demo (make_demo.py).
"""
import argparse
import math
import os
import sys
from pathlib import Path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import pandas as pd
from vmbdiet import (read_humann_table, build_axis_feature_table, logistic_assoc,
                     meta_analyze, or_to_logor)

HERE = Path(__file__).parent


def _binary_outcome(md, outcome, case, control):
    """Return a 0/1 Series: case->1, control->0, everything else NaN (dropped). If no --case was
    given, assume the outcome column is already numeric 0/1."""
    if case is None:
        return pd.to_numeric(md[outcome], errors="coerce")
    s = md[outcome].astype(str).str.strip()
    y = pd.Series(float("nan"), index=md.index, dtype=float)
    y[s == case] = 1.0
    y[s == control] = 0.0
    return y


def per_cohort(label, pathways, meta, outcome, case, control, covars):
    """One cohort: axis capacities -> adjusted per-SD logistic OR of the outcome on each capacity."""
    table = read_humann_table(pathways)
    axis = build_axis_feature_table(table, aggregate="sum")      # samples x {glycogen_starch, scfa,...}
    md = pd.read_csv(meta, sep="\t", index_col=0)
    idx = axis.index.intersection(md.index)
    df = axis.loc[idx].join(md.loc[idx])
    df["_y"] = _binary_outcome(md.loc[idx], outcome, case, control)
    df = df.dropna(subset=["_y"])
    n_case = int((df["_y"] == 1).sum())
    n_ctrl = int((df["_y"] == 0).sum())
    out = {}
    if n_case < 5 or n_ctrl < 5:
        return {"_n_case": n_case, "_n_ctrl": n_ctrl, "_axes": out}   # too thin to model
    use_cov = [c for c in covars if c in df.columns and df[c].notna().any() and df[c].nunique() > 1]
    for cat in axis.columns:
        if df[cat].nunique() < 3:
            continue
        try:
            r = logistic_assoc(df, "_y", cat, covariates=use_cov, standardize=True)
            out[cat] = r
        except Exception as e:  # separation / convergence -> skip this axis for this cohort
            out[cat] = {"note": type(e).__name__}
    return {"_n_case": n_case, "_n_ctrl": n_ctrl, "_covars": use_cov, "_axes": out}


def analyze(cohorts, outcome="disease", case="Type 2 Diabetes Mellitus", control="Healthy",
            covars=("age", "sex", "BMI"), method="DL", verbose=True):
    """cohorts: list of (label, pathways_path, meta_path). Returns a list of per-axis pooled results."""
    covars = list(covars)
    per = {}
    for label, pw, mp in cohorts:
        per[label] = per_cohort(label, pw, mp, outcome, case, control, covars)
    axes = []
    for d in per.values():
        for cat in d["_axes"]:
            if cat not in axes:
                axes.append(cat)

    results = []
    for cat in axes:
        ests, ors, labels = [], [], []
        for label, d in per.items():
            r = d["_axes"].get(cat)
            if not r or "odds_ratio" not in r:
                continue
            ests.append(or_to_logor(r["odds_ratio"], r["ci_low"], r["ci_high"]))
            ors.append((label, r["odds_ratio"], r["p_value"], r["n"]))
            labels.append(label)
        if len(ests) < 2:
            continue
        m = meta_analyze(ests, method=method)
        pooled_or = math.exp(m["beta"])
        same_dir = all(o > 1 for _, o, _, _ in ors) or all(o < 1 for _, o, _, _ in ors)
        replicates = bool(same_dir and m["p_value"] < 0.05)
        results.append(dict(axis=cat, per_cohort=ors, pooled_or=pooled_or,
                            ci_low=math.exp(m["ci"][0]), ci_high=math.exp(m["ci"][1]),
                            p_value=m["p_value"], I2=m["I2"], k=m["k"], replicates=replicates))

    if verbose:
        _print(per, results, outcome, case, method)
    return results


def _print(per, results, outcome, case, method):
    print("=" * 82)
    tot_case = sum(d["_n_case"] for d in per.values())
    tot_ctrl = sum(d["_n_ctrl"] for d in per.values())
    print(f"Gut function <-> {case or outcome}: cross-cohort reproducibility "
          f"({len(per)} cohorts, {tot_case} cases / {tot_ctrl} controls)")
    print("=" * 82)
    for label, d in per.items():
        cov = ",".join(d.get("_covars", [])) or "none"
        print(f"  {label:<12} n={d['_n_case']}+{d['_n_ctrl']}  adj: {cov}")
    if not results:
        print("  (no axis had >=2 modelable cohorts)")
        return
    print(f"\n  {'axis':<16}{'per-cohort OR/SD':<30}{'POOLED OR [95% CI]':<22}{'p':>8}{'I2':>6}  repl")
    for r in results:
        pc = " ".join(f"{o:.2f}" for _, o, _, _ in r["per_cohort"])
        ci = f"{r['ci_low']:.2f}-{r['ci_high']:.2f}"
        flag = "YES" if r["replicates"] else "no"
        print(f"  {r['axis']:<16}{pc:<30}{r['pooled_or']:>5.2f} [{ci:<11}]{r['p_value']:>8.3f}"
              f"{r['I2']:>5.0f}%  {flag}")
    print(f"\n  per-SD adjusted logistic within cohort; pooled by inverse-variance ({method}). "
          f"'repl' = same direction in every cohort AND pooled p<0.05.")
    print("  ASSOCIATION, not mediation (no diet exposure in cMD); relative-abundance capacity over a "
          "small\n  pathway set -- hypothesis-generating. estrobolome needs the gene-family (EC/KO) export.")


def run_demo(verbose=True):
    from make_demo import main as gen
    cohorts = gen()
    # demo metadata already carries a numeric T2D 0/1 column and a disease string; use the string path
    return analyze(cohorts, outcome="disease", case="Type 2 Diabetes Mellitus", control="Healthy",
                   covars=("age", "sex", "BMI"), verbose=verbose)


def main(argv=None):
    p = argparse.ArgumentParser(description="cross-cohort gut-function <-> cardiometabolic reproducibility")
    p.add_argument("--cohort", nargs=3, action="append", metavar=("LABEL", "PATHWAYS", "META"),
                   help="repeatable: a cohort's label, HUMAnN pathway table, and cMD metadata table")
    p.add_argument("--outcome", default="disease", help="metadata column holding the outcome")
    p.add_argument("--case", default=None,
                   help="value counted as a case (=1); OMIT if the outcome column is already numeric 0/1 "
                        "(e.g. the t2d_case column written by export_cmd_cohorts.R)")
    p.add_argument("--control", default="Healthy", help="value counted as a control (=0)")
    p.add_argument("--covars", default="age,sex,BMI", help="comma-separated covariate columns")
    a = p.parse_args(argv)
    covars = tuple(c.strip() for c in a.covars.split(",") if c.strip())
    if not a.cohort:
        run_demo()
        return 0
    analyze([tuple(c) for c in a.cohort], a.outcome, a.case, a.control, covars)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
