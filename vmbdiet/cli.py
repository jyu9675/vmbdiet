"""vmbdiet command-line interface."""
from __future__ import annotations
import argparse
import json
import os
import sys
import pandas as pd

from . import (__version__, assign_cst, lactobacillus_fraction, logistic_scan,
               nutrient_taxon_corr, mannwhitney_by_group, permanova, bray_curtis)


def _read(path):
    return pd.read_csv(path, sep="\t", index_col=0)


def _cmd_cst(a) -> int:
    abund = _read(a.abund)
    if a.valencia:
        from .cst import assign_cst_valencia, load_valencia_centroids
        cst = assign_cst_valencia(abund, load_valencia_centroids())
    else:
        cst = assign_cst(abund, dom_threshold=a.threshold)
    cst.to_csv(a.out, sep="\t")
    vc = cst["cst"].value_counts().reindex(["I", "II", "III", "V", "IV"]).fillna(0).astype(int)
    print("=" * 50)
    print("vmbdiet — community state types")
    print("=" * 50)
    for k, v in vc.items():
        print(f"  CST {k:<3} {v:>4}  ({100*v/len(cst):.0f}%)")
    if "lacto_dominant" in cst.columns:
        print(f"  Lactobacillus-dominant: {int(cst['lacto_dominant'].sum())}/{len(cst)}")
    print(f"  classifier: {'VALENCIA centroids' if a.valencia else 'dominant-taxon'}")
    print(f"[out] {a.out}")
    return 0


def _cmd_associate(a) -> int:
    abund = _read(a.abund)
    meta = _read(a.meta)
    os.makedirs(a.out, exist_ok=True)

    # derive microbiome outcomes from abundance
    cst = assign_cst(abund, dom_threshold=a.threshold)
    df = meta.join(cst[["cst", "lacto_frac", "lacto_dominant"]], how="inner")
    df["cst_IV"] = (df["cst"] == "IV").astype(int)
    df["lacto_depleted"] = (~df["lacto_dominant"].astype(bool)).astype(int)
    cst.to_csv(os.path.join(a.out, "cst.tsv"), sep="\t")

    outcome = a.outcome
    if outcome not in df.columns:
        print(f"ERROR: outcome '{outcome}' not in metadata/derived columns: {list(df.columns)}")
        return 2
    nutrients = [c.strip() for c in a.nutrients.split(",") if c.strip()]
    covars = [c.strip() for c in a.covariates.split(",") if c.strip()] if a.covariates else []

    # 1) adjusted logistic: outcome ~ each nutrient + covariates
    scan = logistic_scan(df, outcome, nutrients, covars)
    scan.to_csv(os.path.join(a.out, "logistic_scan.tsv"), sep="\t", index=False)

    # 2) nutrient ↔ taxon correlations (first nutrient, or all)
    corr_frames = []
    for n in nutrients:
        if n in df.columns:
            c = nutrient_taxon_corr(df[n], abund, method=a.corr_method)
            if not c.empty:
                c.insert(0, "nutrient", n); corr_frames.append(c)
    if corr_frames:
        pd.concat(corr_frames).to_csv(os.path.join(a.out, "nutrient_taxon_corr.tsv"), sep="\t", index=False)

    # 3) PERMANOVA: community (Bray-Curtis) ~ outcome grouping
    perm = None
    try:
        D = bray_curtis(abund.loc[df.index])
        perm = permanova(D, df[outcome].astype(str))
    except Exception as e:
        perm = {"error": str(e)[:80]}

    summary = {"n": int(len(df)), "outcome": outcome, "covariates": covars,
               "cst_distribution": df["cst"].value_counts().to_dict(),
               "permanova_outcome": perm,
               "top_hits": scan.head(5).to_dict(orient="records")}
    with open(os.path.join(a.out, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)

    print("=" * 62)
    print(f"vmbdiet — diet association with '{outcome}'  (n={len(df)})")
    print("=" * 62)
    print(f"  CST: " + ", ".join(f"{k}:{v}" for k, v in sorted(df['cst'].value_counts().items())))
    if isinstance(perm, dict) and "pseudo_F" in perm:
        print(f"  PERMANOVA community~{outcome}: pseudo-F={perm['pseudo_F']:.2f}, p={perm['p_value']:.3f}")
    print(f"\n  adjusted logistic ({outcome} ~ nutrient + {len(covars)} covariates), per-SD OR:")
    print(f"  {'nutrient':<22}{'OR':>7} {'95% CI':>16} {'p':>8} {'q':>8}")
    for _, r in scan.iterrows():
        if pd.isna(r.get("odds_ratio")):
            print(f"  {r['predictor']:<22}  (skipped: {r.get('note','')})"); continue
        ci = f"{r['ci_low']:.2f}-{r['ci_high']:.2f}"
        q = r.get("q_value", float("nan"))
        star = " *" if r["p_value"] < 0.05 else ""
        print(f"  {r['predictor']:<22}{r['odds_ratio']:>7.2f} {ci:>16} {r['p_value']:>8.3f} {q:>8.3f}{star}")
    print(f"\n[out] {a.out}/  (cst.tsv, logistic_scan.tsv, nutrient_taxon_corr.tsv, summary.json)")
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="vmbdiet",
        description="Diet ↔ vaginal microbiome associations (CST + nutrient/food tests).")
    sub = p.add_subparsers(dest="cmd")

    c = sub.add_parser("cst", help="assign community state types from a 16S abundance table")
    c.add_argument("--abund", required=True, help="TSV: samples (rows) × taxa (cols), counts or relabund")
    c.add_argument("--out", default="cst.tsv")
    c.add_argument("--threshold", type=float, default=0.50, help="dominance threshold (default 0.50)")
    c.add_argument("--valencia", action="store_true", help="use packaged VALENCIA centroids (nearest-centroid) instead of dominant-taxon")
    c.set_defaults(func=_cmd_cst)

    a = sub.add_parser("associate", help="test diet/nutrient associations with CST/BV/Lactobacillus")
    a.add_argument("--abund", required=True, help="TSV: samples × taxa abundance")
    a.add_argument("--meta", required=True, help="TSV: samples × (diet/nutrients + covariates [+ outcome])")
    a.add_argument("--outcome", required=True,
                   help="0/1 column to model; derived options: cst_IV, lacto_depleted")
    a.add_argument("--nutrients", required=True, help="comma-separated nutrient/diet columns")
    a.add_argument("--covariates", default="", help="comma-separated covariate columns (age,race,bmi,...)")
    a.add_argument("--corr-method", default="spearman", choices=["spearman", "pearson"])
    a.add_argument("--threshold", type=float, default=0.50)
    a.add_argument("--out", default="vmbdiet_out")
    a.set_defaults(func=_cmd_associate)

    v = sub.add_parser("version")
    v.set_defaults(func=lambda _a: (print(f"vmbdiet {__version__}") or 0))

    args = p.parse_args(argv)
    if not getattr(args, "cmd", None):
        p.print_help(); return 1
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
