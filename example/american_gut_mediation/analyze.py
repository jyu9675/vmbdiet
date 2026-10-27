#!/usr/bin/env python3
"""American Gut diet -> gut -> health MEDIATION prototype (Route B, fully public).

Tests whether diet quality acts on BMI *through* the gut microbiome, on the one fully-open substrate
that carries detailed diet + microbiome + a health outcome with no application (American Gut): the
public VioScreen FFQ (HEI-2010) + 16S alpha diversity + self-reported BMI. It runs the counterfactual
mediation in vmbdiet.mediation:

        diet quality (HEI)  ->  gut alpha diversity (Shannon)  ->  BMI

decomposing diet's total effect on BMI into the part flowing through diversity (natural indirect
effect, NIE) and the rest (natural direct effect, NDE), adjusted for age + sex.

HONEST LIMITS (this is a hypothesis-generating prototype, not established causation):
  * cross-sectional, one sample per participant -> the temporal order diet->diversity->BMI is assumed,
    not shown; the estimate is associational (see vmbdiet.mediation assumptions).
  * BMI and diet are self-reported; the microbiome is 16S, so diversity is the mediator (a measured,
    robust summary) rather than measured metabolic FUNCTION (that needs shotgun / HUMAnN).
  * American Gut is a WEIRD, self-selected citizen-science cohort; the VioScreen FFQ is a subset.

Get the input: `python fetch_americangut.py` (writes agp_mediation_input.tsv from public AGP data).
Then: `python analyze.py --data agp_mediation_input.tsv`  (add --healthy-only to restrict to AGP's
healthy_subset). With no --data it runs the committed synthetic demo (make_demo.py).
"""
import argparse
import os
import sys
from pathlib import Path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import pandas as pd
from vmbdiet import mediate

HERE = Path(__file__).parent


def run(df, n_boot=500, seed=0, label="American Gut"):
    for c in ("hei", "shannon", "bmi", "age"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["hei", "shannon", "bmi", "age", "sex"])
    r = mediate(df, exposure="hei", mediator="shannon", outcome="bmi",
                covariates=["age", "sex"], n_boot=n_boot, seed=seed)

    def ci(k):
        lo, hi = r[k]
        return f"[{lo:+.3f}, {hi:+.3f}]"
    print("=" * 78)
    print(f"{label}: does diet quality (HEI) act on BMI THROUGH gut diversity (Shannon)?")
    print(f"  n={r['n']}  outcome={r['outcome_type']}  effects per +1 SD HEI, adjusted age+sex")
    print("=" * 78)
    print(f"  path a (HEI -> Shannon)     : {r['path_a_x_to_m']:+.4f}")
    print(f"  path b (Shannon -> BMI|HEI) : {r['path_b_m_to_y']:+.4f}")
    print(f"  total effect (HEI -> BMI)   : {r['total_effect']:+.3f}  {ci('total_ci')}")
    print(f"  direct  (NDE, not via div.) : {r['direct_effect_nde']:+.3f}  {ci('nde_ci')}")
    print(f"  indirect(NIE, via diversity): {r['indirect_effect_nie']:+.3f}  {ci('nie_ci')}")
    pm = r["proportion_mediated"]
    print(f"  proportion mediated         : {pm:.1%}  {ci('prop_ci')}"
          if pm == pm else "  proportion mediated         : n/a")
    sig = "excludes 0 -> consistent with mediation" if (r['nie_ci'][0] > 0) == (r['nie_ci'][1] > 0) \
        else "includes 0 -> no evidence of mediation"
    print(f"\n  NIE 95% CI {sig} (under the temporal-order + no-unmeasured-confounding assumptions).")
    print("  Cross-sectional, self-reported BMI, 16S diversity mediator -> hypothesis-generating.")
    return r


def main(argv=None):
    p = argparse.ArgumentParser(description="American Gut diet->gut->BMI mediation prototype")
    p.add_argument("--data", help="TSV from fetch_americangut.py; omit to run the synthetic demo")
    p.add_argument("--healthy-only", action="store_true",
                   help="restrict to AGP's healthy_subset (age 20-69, BMI 18.5-30, no IBD/diabetes/abx)")
    p.add_argument("--n-boot", type=int, default=500)
    a = p.parse_args(argv)
    if a.data:
        df = pd.read_csv(a.data, sep="\t")
        label = "American Gut"
        if a.healthy_only and "healthy_subset" in df.columns:
            df = df[df["healthy_subset"].astype(str).str.lower().isin(["yes", "true", "1"])]
            label += " (healthy subset)"
    else:
        from make_demo import make
        df = make()
        label = "SYNTHETIC demo"
    run(df, n_boot=a.n_boot, label=label)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
