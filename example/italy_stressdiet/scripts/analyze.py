#!/usr/bin/env python
"""
Build per-sample taxon table from vsearch ASV output, assign species-level labels to
Lactobacillus (Tier-2), join to the FFQ nutrient table, and run vmbdiet CST + diet associations.
Reproduces the Italian StressDiet cohort (PRJNA1188525 / fcimb.2025.1582283) analysis.
"""
import os, re, subprocess, json, sys
import numpy as np, pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(BASE, "work")
VS = os.path.join(BASE, "tools", "vsearch-2.32.0-win-x86_64", "bin", "vsearch.exe")
OUT = os.path.join(BASE, "results"); os.makedirs(OUT, exist_ok=True)

LACTO_SPP = ["crispatus", "iners", "gasseri", "jensenii"]  # the CST-defining Lactobacillus

def read_fasta(fp):
    d = {}; n = None; b = []
    for l in open(fp):
        if l.startswith(">"):
            if n: d[n] = "".join(b)
            n = l[1:].split()[0].strip(); b = []
        else: b.append(l.strip())
    if n: d[n] = "".join(b)
    return d

def species_assign():
    """usearch_global ASVs vs curated vaginal species DB; return {ASV: (species,ident)} for top hits."""
    hit = os.path.join(WORK, "asv_species.txt")
    subprocess.run([VS, "--usearch_global", os.path.join(WORK, "zotus.fa"),
                    "--db", os.path.join(BASE, "ref", "species_db.fa"),
                    "--id", "0.95", "--top_hits_only", "--maxaccepts", "8", "--maxrejects", "32",
                    "--strand", "both", "--userout", hit,
                    "--userfields", "query+target+id", "--quiet"], check=True)
    best = {}
    for line in open(hit):
        q, t, pid = line.rstrip("\n").split("\t")
        pid = float(pid); sp = t.split("|")[0]
        if q not in best or pid > best[q][1]:
            best[q] = (sp, pid)
    return best

def genus_from_sintax():
    """parse vsearch --sintax tabbedout -> {ASV: genus} using the cutoff-filtered column."""
    g = {}
    for line in open(os.path.join(WORK, "sintax.txt")):
        parts = line.rstrip("\n").split("\t")
        asv = parts[0]
        taxstr = parts[3] if len(parts) > 3 and parts[3] else (parts[1] if len(parts) > 1 else "")
        m = re.search(r"g:([^,]+)", taxstr)
        fam = re.search(r"f:([^,]+)", taxstr)
        g[asv] = m.group(1) if m else (("f_" + fam.group(1)) if fam else "Unclassified")
    return g

def main():
    otu = pd.read_csv(os.path.join(WORK, "otutab.txt"), sep="\t", index_col=0)
    otu.index.name = "ASV"
    seqs = read_fasta(os.path.join(WORK, "zotus.fa"))
    genus = genus_from_sintax()
    spp = species_assign()

    # label each ASV: species-level for Lactobacillus (>=99% conf), else genus
    labels = {}
    for asv in otu.index:
        sp, pid = spp.get(asv, (None, 0))
        g = genus.get(asv, "Unclassified")
        if sp and pid >= 99.0 and sp.split("_")[0] in ("Lactobacillus", "Limosilactobacillus"):
            labels[asv] = sp                       # e.g. Lactobacillus_crispatus
        elif sp and pid >= 99.0 and g in ("Gardnerella", "Fannyhessea", "Atopobium"):
            labels[asv] = sp
        elif g and g != "Unclassified":
            labels[asv] = g
        elif sp and pid >= 97.0:
            labels[asv] = sp
        else:
            labels[asv] = "Unclassified"
    lab = pd.Series(labels, name="taxon")

    # ASV x sample -> taxon x sample -> sample x taxon
    tab = otu.copy(); tab["taxon"] = tab.index.map(labels)
    taxon_tab = tab.groupby("taxon").sum().T      # samples x taxa (counts)
    taxon_tab.index.name = "sample_alias"
    taxon_tab.to_csv(os.path.join(OUT, "abundance_counts.tsv"), sep="\t")

    # drop Unclassified for CST (keep for reporting), require min depth
    depth = taxon_tab.sum(axis=1)
    keep = depth[depth >= 1000].index
    dropped = sorted(set(taxon_tab.index) - set(keep))
    A = taxon_tab.loc[keep]
    print(f"[taxa] {taxon_tab.shape[1]} taxa; samples passing depth>=1000: {len(keep)}/{len(taxon_tab)}"
          + (f"  (dropped low-depth: {dropped})" if dropped else ""))

    # ---- vmbdiet CST ----
    from vmbdiet import assign_cst, lactobacillus_fraction
    cst = assign_cst(A, dom_threshold=0.50)
    cst.to_csv(os.path.join(OUT, "cst.tsv"), sep="\t")
    vc = cst["cst"].value_counts().reindex(["I","II","III","V","IV"]).fillna(0).astype(int)
    print("[CST]", ", ".join(f"{k}:{v}" for k,v in vc.items()),
          f"| Lactobacillus-dominant {int(cst['lacto_dominant'].sum())}/{len(cst)}")

    # ---- join to nutrients ----
    nut = pd.read_csv(os.path.join(BASE, "nutrients.tsv"), sep="\t")
    nut["sample_alias"] = "StressDiet_" + nut["sample"].astype(int).astype(str)
    nut = nut.set_index("sample_alias")
    df = nut.join(cst[["cst","lacto_frac","lacto_dominant"]], how="inner")
    df["cst_IV"] = (df["cst"] == "IV").astype(int)
    df["lacto_depleted"] = (~df["lacto_dominant"].astype(bool)).astype(int)
    print(f"[join] N with diet+microbiome = {len(df)}  | CST-IV {int(df.cst_IV.sum())}  "
          f"lacto-depleted {int(df.lacto_depleted.sum())}")

    # relative-abundance taxa for correlations (incl Gardnerella)
    rel = A.div(A.sum(axis=1), axis=0)
    rel.to_csv(os.path.join(OUT, "abundance_relabund.tsv"), sep="\t")

    # ---- associations (vmbdiet) ----
    from vmbdiet import logistic_scan, nutrient_taxon_corr, permanova, bray_curtis
    nutrients = ["animal proteins","veg. Proteins","total lipids","animal lipids","saturated FA",
                 "monoins. FA","linoleic acid","other PUFA","starch","carbo","fiber","Alcohol","MEDILITE"]
    nutrients = [n for n in nutrients if n in df.columns]
    covars = ["Energy kcal"]   # energy-adjusted, as in the paper

    results = {}
    for outcome in ["cst_IV", "lacto_depleted"]:
        scan = logistic_scan(df, outcome, nutrients, covars)
        scan.to_csv(os.path.join(OUT, f"logistic_{outcome}.tsv"), sep="\t", index=False)
        results[outcome] = scan
        try:
            D = bray_curtis(rel.loc[df.index]); perm = permanova(D, df[outcome].astype(str))
        except Exception as e:
            perm = {"error": str(e)[:80]}
        results[f"perm_{outcome}"] = perm

    # Gardnerella / key taxa vs alcohol & animal protein (paper: alcohol->Gardnerella, animal protein->CST IV)
    corr_frames = []
    for n in ["Alcohol","animal proteins","MEDILITE","fiber"]:
        if n in df.columns:
            c = nutrient_taxon_corr(df[n], rel.loc[df.index], method="spearman")
            if not c.empty:
                c.insert(0, "nutrient", n); corr_frames.append(c)
    corr = pd.concat(corr_frames) if corr_frames else pd.DataFrame()
    if not corr.empty:
        corr.to_csv(os.path.join(OUT, "nutrient_taxon_corr.tsv"), sep="\t", index=False)

    # ---- print summary ----
    print("\n" + "="*72)
    print("DIET -> VAGINAL MICROBIOME ASSOCIATIONS  (Italian StressDiet cohort, N=%d)" % len(df))
    print("="*72)
    for outcome in ["cst_IV", "lacto_depleted"]:
        perm = results[f"perm_{outcome}"]
        print(f"\n### outcome = {outcome}  (events={int(df[outcome].sum())}/{len(df)})")
        if isinstance(perm, dict) and "pseudo_F" in perm:
            print(f"  PERMANOVA community~{outcome}: pseudo-F={perm['pseudo_F']:.2f} p={perm['p_value']:.3f}")
        sc = results[outcome].copy()
        print(f"  {'nutrient':<20}{'OR/SD':>7} {'95%CI':>14} {'p':>8} {'q':>8}")
        for _, r in sc.iterrows():
            if pd.isna(r.get("odds_ratio")):
                print(f"  {r['predictor']:<20}  (skipped: {r.get('note','')})"); continue
            star = " *" if r["p_value"] < 0.05 else ""
            print(f"  {r['predictor']:<20}{r['odds_ratio']:>7.2f} "
                  f"{r['ci_low']:>6.2f}-{r['ci_high']:<6.2f} {r['p_value']:>8.3f} {r.get('q_value',float('nan')):>8.3f}{star}")
    if not corr.empty:
        print("\n### nutrient -> taxon (Spearman, top |rho| per key nutrient)")
        for n in corr["nutrient"].unique():
            sub = corr[corr["nutrient"]==n].copy()
            sub = sub.reindex(sub["r"].abs().sort_values(ascending=False).index).head(3)
            for _, row in sub.iterrows():
                star = " *" if row["p_value"]<0.05 else ""
                print(f"  {n:<16} {row['taxon']:<24} rho={row['r']:+.2f} p={row['p_value']:.3f}{star}")

    # ---- energy-adjustment: diet QUALITY (per 1000 kcal density) vs CST-IV ----
    from vmbdiet import mannwhitney_by_group
    E = df["Energy kcal"] / 1000.0
    dens_cols = []
    for n in ["animal proteins","veg. Proteins","animal lipids","saturated FA","fiber","carbo","Alcohol"]:
        if n in df.columns:
            df[n + "_dens"] = df[n] / E; dens_cols.append(n + "_dens")
    mw_rows = []
    for n in ["animal proteins","veg. Proteins","animal lipids","saturated FA","fiber","carbo","Alcohol",
              "MEDILITE","Energy kcal"] + dens_cols:
        if n not in df.columns: continue
        r = mannwhitney_by_group(df[n], df["cst_IV"].map({1: "CST_IV", 0: "other"}))
        med_iv = r["median_a"] if r["group_a"] == "CST_IV" else r["median_b"]
        med_ot = r["median_b"] if r["group_a"] == "CST_IV" else r["median_a"]
        mw_rows.append(dict(nutrient=n, median_cstIV=med_iv, median_other=med_ot, p_value=r["p_value"]))
    mw = pd.DataFrame(mw_rows); mw.to_csv(os.path.join(OUT, "nutrient_by_cstIV_mannwhitney.tsv"), sep="\t", index=False)
    dens_scan = logistic_scan(df, "cst_IV", dens_cols, ())
    dens_scan.to_csv(os.path.join(OUT, "logistic_cstIV_density.tsv"), sep="\t", index=False)

    print("\n" + "="*72)
    print("ENERGY-ADJUSTED (diet quality, per 1000 kcal) -> CST-IV  [recovers paper's signal]")
    print("="*72)
    print("  Mann-Whitney on nutrient DENSITY (CST-IV vs other):")
    for _, r in mw[mw.nutrient.str.endswith("_dens")].iterrows():
        print(f"    {r['nutrient']:<22} IV={r['median_cstIV']:6.1f} other={r['median_other']:6.1f} "
              f"p={r['p_value']:.3f}{' *' if r['p_value']<0.05 else ''}")
    print("  logistic cst_IV ~ density (per-SD OR):")
    for _, r in dens_scan.iterrows():
        if pd.isna(r.get("odds_ratio")): continue
        print(f"    {r['predictor']:<22} OR={r['odds_ratio']:.2f} ({r['ci_low']:.2f}-{r['ci_high']:.2f}) "
              f"p={r['p_value']:.3f} q={r.get('q_value',float('nan')):.3f}{' *' if r['p_value']<0.05 else ''}")

    # save summary json
    ap = dens_scan[dens_scan.predictor == "animal proteins_dens"]
    gard = corr[(corr.nutrient == "Alcohol") & (corr.taxon.isin(["Gardnerella","Gardnerella_vaginalis"]))] if not corr.empty else pd.DataFrame()
    summ = {"N": int(len(df)),
            "replication": {
              "animal_protein_density_to_cstIV": (ap.iloc[0][["odds_ratio","ci_low","ci_high","p_value","q_value"]].to_dict() if len(ap) else None),
              "alcohol_to_Gardnerella": (gard.sort_values("p_value").iloc[0][["taxon","r","p_value","q_value"]].to_dict() if len(gard) else None)},
            "cst_distribution": cst["cst"].value_counts().to_dict(),
            "cst_IV": int(df.cst_IV.sum()), "lacto_depleted": int(df.lacto_depleted.sum()),
            "n_taxa": int(A.shape[1]), "n_asvs": int(otu.shape[0]),
            "permanova_cstIV": results["perm_cst_IV"],
            "top_cstIV": results["cst_IV"].head(6).to_dict("records")}
    json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=2, default=str)
    print(f"\n[out] {OUT}")

if __name__ == "__main__":
    main()
