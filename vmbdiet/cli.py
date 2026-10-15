"""vmbdiet command-line interface."""
from __future__ import annotations
import argparse
import json
import os
import sys
import pandas as pd

from . import (__version__, assign_cst, lactobacillus_fraction, logistic_scan,
               nutrient_taxon_corr, mannwhitney_by_group, permanova, bray_curtis,
               maaslin_py, prepare_inputs, write_maaslin3_script, build_formula,
               read_humann_table, detect_level, match_features, build_axis_feature_table,
               CATEGORIES, FUNCTIONAL_FEATURES, mediate)


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


def _cmd_maaslin(a) -> int:
    feat = _read(a.features)
    meta = _read(a.meta)
    os.makedirs(a.out, exist_ok=True)
    predictors = [c.strip() for c in a.predictors.split(",") if c.strip()]
    covars = [c.strip() for c in a.covariates.split(",") if c.strip()] if a.covariates else []

    # 1) Python random-effects scan (runs now, no R needed)
    res = maaslin_py(feat, meta, predictors=predictors, covariates=covars,
                     group=a.group, body_site=a.body_site, min_prevalence=a.min_prevalence,
                     models=tuple(a.models.split(",")), normalize=a.normalize)
    res.to_csv(os.path.join(a.out, "maaslin_py_results.tsv"), sep="\t", index=False)

    # 2) MaAsLin 3 inputs + R script (run on the cluster for the definitive fit)
    prepare_inputs(feat, meta, a.out)
    write_maaslin3_script(os.path.join(a.out, "run_maaslin3.R"), predictors=predictors,
                          covariates=covars, body_site=a.body_site, group=a.group)

    print("=" * 70)
    print(f"vmbdiet — MaAsLin 3 bridge  (features={feat.shape[1]}, samples={feat.shape[0]})")
    print("=" * 70)
    print(f"  formula: {build_formula(predictors, covars, a.body_site, a.group)}")
    print(f"  predictors: {', '.join(predictors)}   covariates: {', '.join(covars) or '(none)'}")
    for model in [m for m in ("abundance", "prevalence") if m in res.get("model", pd.Series()).unique()]:
        sub = res[(res.model == model) & res.p_value.notna()].head(8)
        print(f"\n  {model} model — top associations (per-SD effect):")
        print(f"  {'feature':<24}{'predictor':<18}{'effect':>8}{'p':>9}{'q':>9}")
        for _, r in sub.iterrows():
            star = " *" if r["p_value"] < 0.05 else ""
            print(f"  {str(r['feature'])[:23]:<24}{str(r['predictor'])[:17]:<18}"
                  f"{r['effect']:>8.3f}{r['p_value']:>9.3f}{r.get('q_value', float('nan')):>9.3f}{star}")
    print(f"\n[out] {a.out}/  (maaslin_py_results.tsv, features.tsv, metadata.tsv, run_maaslin3.R)")
    print("      -> run `Rscript run_maaslin3.R` on the cluster for the definitive MaAsLin 3 fit.")
    return 0


def _cmd_functional(a) -> int:
    table = read_humann_table(a.humann, stratified=a.stratified)
    os.makedirs(a.out, exist_ok=True)
    cats = [c.strip() for c in a.categories.split(",") if c.strip()] if a.categories else list(CATEGORIES)
    level = detect_level(table["feature_id"]) if a.level == "auto" else a.level

    matched = match_features(table, categories=cats, level=level)
    matched.to_csv(os.path.join(a.out, "matched_features.tsv"), sep="\t", index=False)
    axis = build_axis_feature_table(table, categories=cats, level=level, aggregate=a.aggregate)
    axis.to_csv(os.path.join(a.out, "axis_feature_table.tsv"), sep="\t")

    print("=" * 66)
    print(f"vmbdiet — functional axis features  (level: {level}, samples: {axis.shape[0]})")
    print("=" * 66)
    n_feat = len(read_humann_table(a.humann, stratified=a.stratified))
    print(f"  parsed {n_feat} community-total features; matched {len(matched)} to the axis dictionary")
    for cat in cats:
        sub = matched[matched.category == cat]
        spec = FUNCTIONAL_FEATURES[cat]
        print(f"\n  [{cat}] {spec['label']} — {len(sub)} feature(s)")
        for feat in sub["feature"].head(6):
            print(f"      {str(feat)[:66]}")
    print(f"\n  axis feature table ({a.aggregate}): {axis.shape[0]} samples x {axis.shape[1]} features -> {', '.join(map(str, axis.columns))}")
    print(f"\n[out] {a.out}/  (axis_feature_table.tsv, matched_features.tsv)")
    print("      -> feed axis_feature_table.tsv to `vmbdiet maaslin --features ...` for Aim-3 association.")
    return 0


def _cmd_mediate(a) -> int:
    df = _read(a.data)
    covars = [c.strip() for c in a.covariates.split(",") if c.strip()] if a.covariates else []
    r = mediate(df, a.exposure, a.mediator, a.outcome, covariates=covars, group=a.group,
                outcome_type=a.outcome_type, contrast=a.contrast, n_boot=a.n_boot, sims=a.sims)
    if a.out:
        os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
        with open(a.out, "w") as fh:
            json.dump(r, fh, indent=2, default=str)

    def line(name, est, ci):
        return f"  {name:<28}{est:>9.4f}   [{ci[0]:>8.4f}, {ci[1]:>8.4f}]"
    print("=" * 70)
    print(f"vmbdiet — mediation:  {a.exposure}  ->  {a.mediator}  ->  {a.outcome}")
    print("=" * 70)
    print(f"  n={r['n']}  outcome={r['outcome_type']}  scale={r['scale']}  contrast={r['contrast']}")
    print(f"  bootstrap: {r['n_boot_ok']} ok" + (f", clustered by {r['clustered_by']}" if r['clustered_by'] else ", unclustered"))
    print(f"\n  {'effect':<28}{'estimate':>9}   {'95% CI':>20}")
    print(line("total effect", r["total_effect"], r["total_ci"]))
    print(line("direct (NDE)", r["direct_effect_nde"], r["nde_ci"]))
    print(line("indirect / through M (NIE)", r["indirect_effect_nie"], r["nie_ci"]))
    print(line("proportion mediated", r["proportion_mediated"], r["prop_ci"]))
    sig = r["nie_ci"][0] * r["nie_ci"][1] > 0
    print(f"\n  paths: a (X->M)={r['path_a_x_to_m']:.3f}  b (M->Y)={r['path_b_m_to_y']:.3f}  c' (direct)={r['path_cprime_direct']:.3f}")
    print(f"  -> indirect effect {'EXCLUDES' if sig else 'includes'} 0: "
          + ("consistent with mediation through the gut" if sig else "no evidence of mediation") + " (under the stated assumptions)")
    print(f"\n  NOTE: {r['assumptions']}")
    if a.out:
        print(f"[out] {a.out}")
    return 0


def main(argv=None) -> int:
    try:  # make non-ASCII (β, µ, em dash) safe on a cp1252 Windows console
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
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

    m = sub.add_parser("maaslin",
        help="repeated-measures diet<->feature associations (random effect per participant) + MaAsLin 3 R script")
    m.add_argument("--features", required=True, help="TSV: samples x features (species relabund OR HUMAnN3 pathway/gene abundance)")
    m.add_argument("--meta", required=True, help="TSV: samples x (diet predictors + covariates + participant + body_site)")
    m.add_argument("--predictors", required=True, help="comma-separated diet columns to test")
    m.add_argument("--covariates", default="", help="comma-separated covariates (age,race,bmi,antibiotics,...)")
    m.add_argument("--group", default="participant", help="random-intercept column (default: participant)")
    m.add_argument("--body-site", default=None, help="body-site column to adjust for (e.g. body_site)")
    m.add_argument("--min-prevalence", type=float, default=0.1, help="drop features present in < this fraction (default 0.1)")
    m.add_argument("--normalize", default="tss", choices=["tss", "none"],
                   help="tss for a full compositional profile (default); none for pre-aggregated absolute scores (e.g. axis capacity sums)")
    m.add_argument("--models", default="abundance,prevalence", help="which models to fit")
    m.add_argument("--out", default="vmbdiet_maaslin")
    m.set_defaults(func=_cmd_maaslin)

    fn = sub.add_parser("functional",
        help="extract gut-vagina-axis functional features (estrobolome/glycogen/SCFA) from a HUMAnN 3 table")
    fn.add_argument("--humann", required=True, help="HUMAnN 3 gene-family / pathway / EC / KO table (single-sample or merged)")
    fn.add_argument("--categories", default="", help="comma-separated subset of: " + ",".join(CATEGORIES) + " (default all)")
    fn.add_argument("--level", default="auto", choices=["auto", "pathway", "ec", "ko", "uniref"])
    fn.add_argument("--aggregate", default="sum", choices=["sum", "none"], help="sum per category, or keep per feature")
    fn.add_argument("--stratified", default="total", choices=["total", "species", "keep"])
    fn.add_argument("--out", default="vmbdiet_functional")
    fn.set_defaults(func=_cmd_functional)

    md = sub.add_parser("mediate",
        help="does diet act on the vaginal outcome THROUGH a gut mediator? (natural direct/indirect effects)")
    md.add_argument("--data", required=True, help="TSV: samples/participants x (exposure + mediator + outcome + covariates [+ group])")
    md.add_argument("--exposure", required=True, help="diet column (X)")
    md.add_argument("--mediator", required=True, help="gut function column (M), e.g. an axis capacity score")
    md.add_argument("--outcome", required=True, help="vaginal outcome column (Y): binary (cst_IV) or continuous (lacto_frac)")
    md.add_argument("--covariates", default="", help="comma-separated confounders to adjust for")
    md.add_argument("--group", default=None, help="participant column -> cluster bootstrap resamples participants")
    md.add_argument("--outcome-type", default="auto", choices=["auto", "binary", "continuous"])
    md.add_argument("--contrast", type=float, default=1.0, help="exposure change in SD units (default +1 SD)")
    md.add_argument("--n-boot", type=int, default=500)
    md.add_argument("--sims", type=int, default=100, help="Monte-Carlo sims for the binary-outcome estimator")
    md.add_argument("--out", default=None, help="optional JSON output path")
    md.set_defaults(func=_cmd_mediate)

    v = sub.add_parser("version")
    v.set_defaults(func=lambda _a: (print(f"vmbdiet {__version__}") or 0))

    args = p.parse_args(argv)
    if not getattr(args, "cmd", None):
        p.print_help(); return 1
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
