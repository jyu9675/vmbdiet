#!/usr/bin/env python3
"""Within-person daily diet -> gut FUNCTIONAL-CAPACITY association (Johnson 2019 "Daily Sampling").

The open question the American Gut mediation left behind was whether gut *function* — not the coarse
alpha-diversity summary — tracks diet. This is the one fully-public dataset that can answer it: daily
ASA24 food records paired with daily SHOTGUN metagenomes, so the readout is measured KEGG functional
capacity. It is a within-person, repeated-measures design (not a BMI-outcome cohort): BMI and blood
markers are per-person (n~34), so this tests DIET -> FUNCTION dynamics, not diet -> BMI.

Method (reproducibility-first, mirroring the cross-cohort harness, but the "unit" is a PERSON):
  * for each person, correlate a daily diet variable with a daily axis capacity across their ~17 days
    (within-person Pearson r; the person is their own control, so stable between-person differences in
    body size, transit, enterotype cannot confound);
  * pool the per-person correlations across people by Fisher-z inverse-variance meta-analysis
    (vmbdiet.meta_analyze, DerSimonian-Laird random effects) -> a pooled within-person effect, its CI,
    and I^2. Here HIGH I^2 is itself the finding the original paper reported: diet-microbiome coupling
    is PERSONALIZED, so a population-average effect can be weak even when individuals show strong ones.

Axes (vmbdiet.functional, matched on KEGG enzyme EC numbers):
  * estrobolome      — beta-glucuronidase (EC 3.2.1.31). Uniquely available here at gene level; the
                       cMD gene_families export could not deliver it cheaply (6M ID-only UniRef90 rows).
  * scfa             — short-chain fatty-acid / fermentation enzymes (the fiber mechanism).
  * glycogen_starch  — amylase / glucosidase capacity.

Get the input: `python fetch_dietstudy.py` (writes data/diet_function_input.tsv from the public repo).
Then: `python analyze.py --data data/diet_function_input.tsv`  (--diet-var fiber_dens|fiber|dq_score,
--lag K shifts diet K study-days earlier within person). With no --data it runs the synthetic demo.

HONEST LIMITS: n~34 healthy, WEIRD adults; capacity is EC-level RELATIVE abundance (compositional);
observational (no randomization of diet); the diet record and stool are paired by study-day, so --lag
probes, but cannot prove, directionality. Hypothesis-generating.
"""
import argparse
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import pandas as pd
from scipy import stats
from vmbdiet import meta_analyze, CATEGORIES

HERE = Path(__file__).parent
AXES = list(CATEGORIES)                                   # estrobolome, glycogen_starch, scfa


def _per_person_r(sub, diet_var, axis, lag):
    """Within-person Pearson r between a daily diet variable and a daily axis capacity."""
    s = sub.sort_values("day") if "day" in sub.columns and sub["day"].notna().any() else sub
    x = pd.to_numeric(s[diet_var], errors="coerce")
    y = pd.to_numeric(s[axis], errors="coerce")
    if lag:                                               # diet from `lag` study-days earlier
        x = x.shift(lag)
    d = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(d) < 3 or d["x"].nunique() < 3 or d["y"].nunique() < 3:
        return None
    r, _ = stats.pearsonr(d["x"], d["y"])
    return float(r), int(len(d))


def within_person(df, diet_var, axis, min_days=8, lag=0, method="DL"):
    """Per-person within-person r -> Fisher-z random-effects pool across people."""
    per = []
    for person, sub in df.groupby("person"):
        got = _per_person_r(sub, diet_var, axis, lag)
        if got is None:
            continue
        r, n = got
        if n < min_days:
            continue
        r = max(min(r, 0.999), -0.999)
        z = math.atanh(r)
        se = 1.0 / math.sqrt(n - 3)                       # Fisher-z variance = 1/(n-3)
        per.append((person, r, z, se, n))
    if len(per) < 2:
        return None
    m = meta_analyze([(z, se) for _, _, z, se, _ in per], method=method)
    pooled_r = math.tanh(m["beta"])
    ci = (math.tanh(m["ci"][0]), math.tanh(m["ci"][1]))
    n_pos = sum(1 for _, r, _, _, _ in per if r > 0)
    k = len(per)
    same = max(n_pos, k - n_pos)
    consistent = bool(m["p_value"] < 0.05 and same / k >= 0.6
                      and ((pooled_r > 0) == (n_pos >= k - n_pos)))
    return dict(axis=axis, diet_var=diet_var, k=k, median_days=int(np.median([n for *_, n in per])),
                pooled_r=pooled_r, ci_low=ci[0], ci_high=ci[1], p_value=m["p_value"], I2=m["I2"],
                frac_same_sign=same / k, n_pos=n_pos, consistent=consistent, per_person=per)


def analyze(df, diet_var="fiber_dens", min_days=8, lag=0, method="DL", verbose=True):
    results = [r for r in (within_person(df, diet_var, ax, min_days, lag, method)
                           for ax in AXES if ax in df.columns) if r]
    if verbose:
        _print(df, results, diet_var, lag, method)
    return results


def _print(df, results, diet_var, lag, method):
    print("=" * 86)
    print(f"Within-person daily diet ({diet_var}) -> gut functional capacity  "
          f"(Johnson 2019; lag={lag} study-days)")
    npeople = df["person"].nunique() if "person" in df.columns else 0
    print(f"  {npeople} people, {len(df)} person-days; unit = PERSON (each is their own control)")
    print("=" * 86)
    if not results:
        print("  (no axis had >=2 modelable people)")
        return
    print(f"  {'axis':<16}{'k':>3}{'md':>4}  {'pooled within-person r [95% CI]':<30}{'p':>8}{'I2':>6}"
          f"{'%same':>7}  consistent")
    for r in results:
        ci = f"{r['ci_low']:+.2f},{r['ci_high']:+.2f}"
        flag = "YES" if r["consistent"] else "no"
        print(f"  {r['axis']:<16}{r['k']:>3}{r['median_days']:>4}  "
              f"{r['pooled_r']:+.3f} [{ci:<13}]{r['p_value']:>8.3f}{r['I2']:>5.0f}%"
              f"{r['frac_same_sign']*100:>6.0f}%  {flag}")
    print(f"\n  per-person Pearson r pooled by Fisher-z inverse-variance ({method}); 'consistent' = "
          f"pooled p<0.05\n  AND >=60% of people share the pooled sign. HIGH I^2 = personalized "
          "coupling (the dataset's own finding).")
    print("  estrobolome = beta-glucuronidase (EC 3.2.1.31), the gene-level estrobolome read cMD could "
          "not give.\n  Within-person, healthy n~34, EC-level relative-abundance capacity, "
          "observational -> hypothesis-generating.")


def run_demo(verbose=True):
    sys.path.insert(0, str(HERE))
    import make_demo
    df = make_demo.make()
    return analyze(df, diet_var="fiber_dens", min_days=8, lag=0, verbose=verbose)


def main(argv=None):
    p = argparse.ArgumentParser(description="within-person daily diet -> gut functional capacity")
    p.add_argument("--data", help="TSV from fetch_dietstudy.py; omit to run the synthetic demo")
    p.add_argument("--diet-var", default="fiber_dens",
                   help="diet column: fiber_dens (default), fiber, dq_score, kcal, ...")
    p.add_argument("--min-days", type=int, default=8, help="min paired days to include a person")
    p.add_argument("--lag", type=int, default=0, help="shift diet K study-days earlier within person")
    a = p.parse_args(argv)
    if a.data:
        df = pd.read_csv(a.data, sep="\t")
        analyze(df, diet_var=a.diet_var, min_days=a.min_days, lag=a.lag)
    else:
        run_demo()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
