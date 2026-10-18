#!/usr/bin/env python3
"""Aim-3 functional contrast on curatedMetagenomicData: does gut functional CAPACITY differ by
westernization (cMD's diet/lifestyle proxy)?

Reads a HUMAnN 3 pathway table + a cMD metadata table (carrying `westernized`), extracts the gut-vagina
axis capacity scores (glycogen/starch, SCFA — estrobolome needs the gene-family export) with
vmbdiet.functional, and compares each capacity between Westernized and Non-westernized samples
(Mann-Whitney). Expectation: non-westernized (high-fiber) -> higher fiber-fermentation capacity.

Get the real inputs with ONE R step (no further R needed afterwards):
  library(curatedMetagenomicData); # from vmbdiet.functional.write_cmd_export_R()
  se <- curatedMetagenomicData("2021-10-14.<Study>.pathway_abundance", dryrun=FALSE, counts=FALSE)[[1]]
  write.table(as.data.frame(assay(se)),   "cmd_<Study>_pathway_abundance.tsv", sep="\t", col.names=NA)
  write.table(as.data.frame(colData(se)), "cmd_<Study>_metadata.tsv",          sep="\t", col.names=NA)
Then: python analyze.py --pathways cmd_<Study>_pathway_abundance.tsv --meta cmd_<Study>_metadata.tsv

With no arguments it runs the committed synthetic demo (make_demo.py), so the harness is testable offline.
"""
import argparse
import os
import sys
from pathlib import Path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import pandas as pd
from vmbdiet import build_axis_feature_table, read_humann_table, mannwhitney_by_group

HERE = Path(__file__).parent


def run(pathways, meta, group_col="westernized"):
    table = read_humann_table(pathways)                       # HUMAnN pathway table -> community totals
    axis = build_axis_feature_table(table, aggregate="sum")   # samples x {glycogen_starch, scfa, ...}
    md = pd.read_csv(meta, sep="\t", index_col=0)
    idx = axis.index.intersection(md.index)
    axis, md = axis.loc[idx], md.loc[idx]
    groups = md[group_col].astype(str)
    levels = sorted(groups.dropna().unique())
    print("=" * 70)
    print(f"Aim-3 functional contrast by {group_col}  (n={len(idx)}: " +
          ", ".join(f"{g}={int((groups == g).sum())}" for g in levels) + ")")
    print("=" * 70)
    if len(levels) != 2:
        print(f"  need exactly 2 {group_col} groups, got {levels}"); return 1
    print(f"  median capacity by group (Mann-Whitney):")
    print(f"  {'axis':<17}{levels[0][:15]:>16}{levels[1][:15]:>16}{'p(MWU)':>9}")
    for cat in axis.columns:
        r = mannwhitney_by_group(axis[cat], groups)
        star = " *" if r["p_value"] < 0.05 else ""
        print(f"  {cat:<17}{r['median_a']:>16.1f}{r['median_b']:>16.1f}{r['p_value']:>9.3f}{star}")
    print(f"\n  groups: {levels[0]} vs {levels[1]}.  (estrobolome / beta-glucuronidase needs the "
          f"gene-family export - EC/KO - not the pathway table.)")
    return 0


def main(argv=None):
    p = argparse.ArgumentParser(description="cMD westernized-contrast Aim-3 functional analysis")
    p.add_argument("--pathways", help="HUMAnN 3 pathway table (cMD export); omit to use the synthetic demo")
    p.add_argument("--meta", help="cMD metadata table with a `westernized` column")
    p.add_argument("--group", default="westernized")
    a = p.parse_args(argv)
    if not a.pathways:
        from make_demo import main as gen
        gen()
        a.pathways = HERE / "data" / "cmd_demo_pathway_abundance.tsv"
        a.meta = HERE / "data" / "cmd_demo_metadata.tsv"
    return run(a.pathways, a.meta, a.group)


if __name__ == "__main__":
    raise SystemExit(main())
