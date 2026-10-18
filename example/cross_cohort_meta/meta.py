#!/usr/bin/env python3
"""Cross-cohort meta-analysis of diet -> vaginal CST-IV, across the two public cohorts in this repo.

Reads each cohort's already-computed per-SD, energy-adjusted logistic results (committed under
example/*/results/), maps the shared diet features, and pools them with vmbdiet.meta_analyze. Reproduces
the headline finding: animal-protein density -> CST-IV, pooled OR ~1.54 (p~0.017), consistent (I^2=0%).

Pooling is done at the derived CST-IV outcome with standardized per-SD exposures (not raw taxa), because
the two cohorts use different 16S regions. Run: python example/cross_cohort_meta/meta.py
"""
import csv
import os
import sys
import math
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from vmbdiet import meta_analyze, or_to_logor

HERE = os.path.dirname(os.path.abspath(__file__))
EX = os.path.dirname(HERE)
ITALY = os.path.join(EX, "italy_stressdiet", "results", "logistic_cstIV_density.tsv")
MELODY = os.path.join(EX, "melody_ibd", "results", "logistic_cst_IV_density.tsv")

# shared diet features: (label, italy predictor name, melody predictor name)
PAIRS = [("animal protein", "animal proteins_dens", "animal_protein_dens"),
         ("vegetable protein", "veg. Proteins_dens", "vegetable_protein_dens"),
         ("carbohydrate", "carbo_dens", "carbohydrate_dens"),
         ("fiber", "fiber_dens", "fiber_dens"),
         ("alcohol", "Alcohol_dens", "alcohol_dens"),
         ("saturated fat", "saturated FA_dens", "saturated_fat_dens")]


def load(path):
    out = {}
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[r["predictor"]] = {k: float(r[k]) for k in ("odds_ratio", "ci_low", "ci_high", "n", "p_value")}
    return out


def main():
    it, me = load(ITALY), load(MELODY)
    n_it = int(it["animal proteins_dens"]["n"]); n_me = int(me["animal_protein_dens"]["n"])
    rows = []
    for lab, ik, mk in PAIRS:
        est = [or_to_logor(it[ik]["odds_ratio"], it[ik]["ci_low"], it[ik]["ci_high"]),
               or_to_logor(me[mk]["odds_ratio"], me[mk]["ci_low"], me[mk]["ci_high"])]
        m = meta_analyze(est, method="fixed")
        rows.append(dict(feature=lab, italy_or=it[ik]["odds_ratio"], melody_or=me[mk]["odds_ratio"],
                         pooled_or=math.exp(m["beta"]), ci_low=math.exp(m["ci"][0]),
                         ci_high=math.exp(m["ci"][1]), p_value=m["p_value"], I2=m["I2"]))

    out = os.path.join(HERE, "meta_results.tsv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print("=" * 74)
    print(f"Diet -> vaginal CST-IV, pooled across 2 public cohorts (n={n_it}+{n_me}={n_it+n_me})")
    print("=" * 74)
    print(f"  {'feature':<18}{'Italy':>7}{'MELODY':>8}   {'POOLED OR [95% CI]':<22}{'p':>8}{'I2':>6}")
    for r in rows:
        star = " *" if r["p_value"] < 0.05 else ""
        ci = f"{r['ci_low']:.2f}-{r['ci_high']:.2f}"
        print(f"  {r['feature']:<18}{r['italy_or']:>7.2f}{r['melody_or']:>8.2f}   "
              f"{r['pooled_or']:>5.2f} [{ci:<11}]   {r['p_value']:>7.3f}{r['I2']:>5.0f}%{star}")
    print(f"\n[out] {out}  |  fixed-effect inverse-variance; per-SD, energy-adjusted")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
