#!/usr/bin/env python
"""MELODY (PRJNA915128, vaginal) — CST from raw 16S + diet associations via vmbdiet.
Diet (HEI-2015 + nutrients + covariates) comes from the BioSample metadata (diet.tsv).
Primary cross-check: HEI-2015 diet quality -> CST (the Hawaii metric)."""
import os, re, subprocess, json
import numpy as np, pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(BASE, "work")
# Portable: vsearch from PATH (or $VSEARCH); curated species DB shipped in the Italian example
# (or $VMB_SPECIES_DB).
VS = os.environ.get("VSEARCH", "vsearch")
SPP_DB = os.environ.get("VMB_SPECIES_DB",
                        os.path.join(BASE, "..", "..", "italy_stressdiet", "data", "vaginal_species_db.fa"))
OUT = os.path.join(BASE, "results"); os.makedirs(OUT, exist_ok=True)

def read_fasta(fp):
    d={};n=None;b=[]
    for l in open(fp):
        if l.startswith(">"):
            if n:d[n]="".join(b)
            n=l[1:].split()[0].strip();b=[]
        else:b.append(l.strip())
    if n:d[n]="".join(b)
    return d

def species_assign():
    hit=os.path.join(WORK,"asv_species.txt")
    subprocess.run([VS,"--usearch_global",os.path.join(WORK,"zotus.fa"),"--db",SPP_DB,
        "--id","0.95","--top_hits_only","--maxaccepts","8","--maxrejects","32","--strand","both",
        "--userout",hit,"--userfields","query+target+id","--quiet"],check=True)
    best={}
    for line in open(hit):
        q,t,pid=line.rstrip("\n").split("\t");pid=float(pid);sp=t.split("|")[0]
        if q not in best or pid>best[q][1]: best[q]=(sp,pid)
    return best

def genus_from_sintax():
    g={}
    for line in open(os.path.join(WORK,"sintax.txt")):
        p=line.rstrip("\n").split("\t");asv=p[0]
        tx=p[3] if len(p)>3 and p[3] else (p[1] if len(p)>1 else "")
        m=re.search(r"g:([^,]+)",tx);f=re.search(r"f:([^,]+)",tx)
        g[asv]=m.group(1) if m else (("f_"+f.group(1)) if f else "Unclassified")
    return g

def main():
    otu=pd.read_csv(os.path.join(WORK,"otutab.txt"),sep="\t",index_col=0); otu.index.name="ASV"
    genus=genus_from_sintax(); spp=species_assign()
    labels={}
    for asv in otu.index:
        sp,pid=spp.get(asv,(None,0)); g=genus.get(asv,"Unclassified")
        if sp and pid>=99 and sp.split("_")[0] in ("Lactobacillus","Limosilactobacillus"): labels[asv]=sp
        elif sp and pid>=99 and g in ("Gardnerella","Fannyhessea","Atopobium"): labels[asv]=sp
        elif g and g!="Unclassified": labels[asv]=g
        elif sp and pid>=97: labels[asv]=sp
        else: labels[asv]="Unclassified"
    tab=otu.copy(); tab["taxon"]=tab.index.map(labels)
    taxon_tab=tab.groupby("taxon").sum().T; taxon_tab.index.name="sample_alias"
    taxon_tab.to_csv(os.path.join(OUT,"abundance_counts.tsv"),sep="\t")
    depth=taxon_tab.sum(axis=1); A=taxon_tab.loc[depth[depth>=1000].index]
    print(f"[taxa] {taxon_tab.shape[1]} taxa; samples depth>=1000: {A.shape[0]}/{taxon_tab.shape[0]}")

    from vmbdiet import assign_cst, logistic_scan, nutrient_taxon_corr, permanova, bray_curtis, mannwhitney_by_group
    cst=assign_cst(A,dom_threshold=0.50); cst.to_csv(os.path.join(OUT,"cst.tsv"),sep="\t")
    vc=cst["cst"].value_counts().reindex(["I","II","III","V","IV"]).fillna(0).astype(int)
    print("[CST]",", ".join(f"{k}:{v}" for k,v in vc.items()),
          f"| Lacto-dominant {int(cst['lacto_dominant'].sum())}/{len(cst)}")

    diet=pd.read_csv(os.path.join(BASE,"diet.tsv"),sep="\t").set_index("sample_alias")
    rel=A.div(A.sum(axis=1),axis=0); rel.to_csv(os.path.join(OUT,"abundance_relabund.tsv"),sep="\t")
    df=diet.join(cst[["cst","lacto_frac","lacto_dominant"]],how="inner")
    df["cst_IV"]=(df["cst"]=="IV").astype(int)
    df["lacto_depleted"]=(~df["lacto_dominant"].astype(bool)).astype(int)
    df=df.dropna(subset=["HEI2015","energy_kcal"])  # drop the 1 sample lacking diet
    print(f"[join] N with diet+microbiome = {len(df)} | CST-IV {int(df.cst_IV.sum())} | IBD {df['ibd'].value_counts().to_dict()}")

    # nutrient densities (per 1000 kcal) + diet-quality score
    E=df["energy_kcal"]/1000.0
    for n in ["animal_protein","vegetable_protein","fiber","saturated_fat","carbohydrate","alcohol"]:
        df[n+"_dens"]=df[n]/E
    preds_raw=["HEI2015","AHEI","animal_protein","vegetable_protein","fiber","saturated_fat","carbohydrate","alcohol"]
    preds_dens=[n+"_dens" for n in ["animal_protein","vegetable_protein","fiber","saturated_fat","carbohydrate","alcohol"]]

    results={}
    for outcome in ["cst_IV","lacto_depleted"]:
        sc=logistic_scan(df,outcome,preds_raw,["energy_kcal"]); sc.to_csv(os.path.join(OUT,f"logistic_{outcome}.tsv"),sep="\t",index=False)
        scd=logistic_scan(df,outcome,preds_dens,()); scd.to_csv(os.path.join(OUT,f"logistic_{outcome}_density.tsv"),sep="\t",index=False)
        results[outcome]=sc; results[outcome+"_dens"]=scd
        try: perm=permanova(bray_curtis(rel.loc[df.index]),df[outcome].astype(str))
        except Exception as e: perm={"error":str(e)[:80]}
        results["perm_"+outcome]=perm

    # HEI-2015 and animal protein by CST-IV (Mann-Whitney), + alcohol->taxon
    mw={}
    for n in ["HEI2015","AHEI","animal_protein_dens","vegetable_protein_dens","fiber_dens","alcohol"]:
        r=mannwhitney_by_group(df[n],df["cst_IV"].map({1:"IV",0:"other"}))
        mw[n]=(r["median_a"] if r["group_a"]=="IV" else r["median_b"], r["median_b"] if r["group_a"]=="IV" else r["median_a"], r["p_value"])
    corr_frames=[]
    for n in ["alcohol","animal_protein","HEI2015"]:
        c=nutrient_taxon_corr(df[n],rel.loc[df.index],method="spearman")
        if not c.empty: c.insert(0,"nutrient",n); corr_frames.append(c)
    corr=pd.concat(corr_frames) if corr_frames else pd.DataFrame()
    if not corr.empty: corr.to_csv(os.path.join(OUT,"nutrient_taxon_corr.tsv"),sep="\t",index=False)

    print("\n"+"="*70); print(f"MELODY diet -> vaginal microbiome (N={len(df)})"); print("="*70)
    for outcome in ["cst_IV","lacto_depleted"]:
        perm=results["perm_"+outcome]
        print(f"\n### {outcome} (events {int(df[outcome].sum())}/{len(df)})")
        if isinstance(perm,dict) and "pseudo_F" in perm: print(f"  PERMANOVA: pseudo-F={perm['pseudo_F']:.2f} p={perm['p_value']:.3f}")
        print("  energy-adjusted logistic (per-SD OR):")
        for _,r in results[outcome].iterrows():
            if pd.isna(r.get("odds_ratio")): print(f"    {r['predictor']:<20} skip: {r.get('note','')}"); continue
            print(f"    {r['predictor']:<20} OR={r['odds_ratio']:.2f} ({r['ci_low']:.2f}-{r['ci_high']:.2f}) p={r['p_value']:.3f} q={r.get('q_value',np.nan):.3f}{' *' if r['p_value']<0.05 else ''}")
    print("\n### Mann-Whitney by CST-IV (IV vs other):")
    for n,(a,b,p) in mw.items(): print(f"    {n:<22} IV={a:7.1f} other={b:7.1f} p={p:.3f}{' *' if p<0.05 else ''}")
    if not corr.empty:
        print("\n### alcohol -> taxon (Spearman, prevalent):")
        prev=(rel>0).mean()
        a=corr[(corr.nutrient=='alcohol')&(corr.taxon.isin(prev[prev>=0.15].index))].sort_values('p_value').head(5)
        for _,r in a.iterrows(): print(f"    {r['taxon']:<24} rho={r['r']:+.2f} p={r['p_value']:.3f} q={r['q_value']:.3f}{' *' if r['p_value']<0.05 else ''}")

    hei=results["cst_IV"][results["cst_IV"].predictor=="HEI2015"]
    ap=results["cst_IV_dens"][results["cst_IV_dens"].predictor=="animal_protein_dens"]
    json.dump({"N":int(len(df)),"cst_distribution":cst["cst"].value_counts().to_dict(),
        "cst_IV":int(df.cst_IV.sum()),"ibd":df["ibd"].value_counts().to_dict(),
        "HEI2015_to_cstIV":(hei.iloc[0][["odds_ratio","ci_low","ci_high","p_value"]].to_dict() if len(hei) else None),
        "animal_protein_density_to_cstIV":(ap.iloc[0][["odds_ratio","ci_low","ci_high","p_value"]].to_dict() if len(ap) else None),
        "n_asvs":int(otu.shape[0]),"n_taxa":int(A.shape[1])},
        open(os.path.join(OUT,"summary.json"),"w"),indent=2,default=str)
    print(f"\n[out] {OUT}")

if __name__=="__main__": main()
