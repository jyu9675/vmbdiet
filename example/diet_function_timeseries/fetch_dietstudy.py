#!/usr/bin/env python3
"""Build the within-person diet -> gut-FUNCTION input from the public Johnson 2019 "Daily Sampling"
study (Knights lab), the one fully-open dataset that pairs DAILY food records with DAILY SHOTGUN
metagenomes — so the mediator is measured functional CAPACITY (KEGG enzymes), not 16S alpha diversity.

Source repo (fully public, processed tables, no SRA reprocessing):
    https://github.com/knights-lab/dietstudy_analyses
    Johnson AJ et al. "Daily sampling reveals personalized diet-microbiome associations in humans."
    Cell Host & Microbe 25, 789-802 (2019).  (34 healthy adults x 17 consecutive days.)

What it does:
  1. KEGG *enzyme* table (EC-level, features x samples) -> the vmbdiet functional axes
     (estrobolome EC 3.2.1.31, glycogen_starch, scfa) via the EXISTING `ecs` dictionaries in
     vmbdiet.functional — the enzyme table is EC-keyed, so the MetaCyc-pathway mismatch that blocks
     cMD does not arise here; we only parse the bare EC id out of each "EC  name" row label.
  2. ASA24 daily nutrition totals -> interpretable diet variables (dietary fiber, fiber density per
     1000 kcal, and a transparent plant/fiber-forward diet-quality score — NOT HEI-2015).
  3. Joins function and diet on the shared per-study-day key `MCT.f.####` (the function columns carry
     a `.S###.R1.001` sequencing suffix we strip; rerun duplicates are averaged).

Output: a tidy long TSV, one row per (person, study-day): person, day, diet vars, and the three axis
capacities. Feed it to analyze.py. The real data is NOT redistributed (gitignored) — reproduce with:
    git clone https://github.com/knights-lab/dietstudy_analyses example/diet_function_timeseries/data/dietstudy_analyses
    python fetch_dietstudy.py                     # uses data/dietstudy_analyses by default
    python fetch_dietstudy.py --repo /path/to/dietstudy_analyses --out diet_function_input.tsv
"""
import argparse
import os
import re
import sys
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
import numpy as np
import pandas as pd
from vmbdiet import build_axis_feature_table

HERE = Path(__file__).parent
_EC = re.compile(r"^(\d+\.\d+\.\d+\.\d+)")
_SAMPLE = re.compile(r"^(MCT\.f\.\d+)")

# plant/fiber-forward diet-quality score (transparent, density-adjusted; NOT the HEI-2015 algorithm).
# "more is better" vs "less is better" USDA food-group / nutrient columns present in nutrition_totals.
_DQ_GOOD = ["FIBE", "G_WHOLE", "V_TOTAL", "F_TOTAL", "PF_LEGUMES", "PF_SEAFD_HI"]
_DQ_BAD = ["ADD_SUGARS", "G_REFINED", "SFAT", "SODI"]


def _ec_axes(enzyme_path: Path) -> pd.DataFrame:
    """KEGG enzyme table (features x samples, row label 'EC  name') -> samples x axis capacities.
    Parses the bare EC id so vmbdiet.functional's EC dictionaries match; averages rerun duplicates."""
    raw = pd.read_csv(enzyme_path, sep="\t", header=0, index_col=0)
    ec = raw.index.to_series().str.extract(_EC, expand=False)
    raw = raw[ec.notna().to_numpy()].copy()
    raw.index = ec.dropna().to_numpy()               # index is now the bare EC number
    raw.columns = [(_SAMPLE.match(c).group(1) if _SAMPLE.match(c) else c) for c in raw.columns]
    axes = build_axis_feature_table(raw, aggregate="sum")        # samples x {estrobolome, ...}
    axes = axes.apply(pd.to_numeric, errors="coerce")
    return axes.groupby(level=0).mean()                          # collapse sequencing reruns


def _diet_vars(diet_path: Path) -> pd.DataFrame:
    """ASA24 nutrition_totals -> per (MCT.f.####) diet variables, averaging any replicate records."""
    d = pd.read_csv(diet_path, sep="\t", header=0)
    key = "X.SampleID" if "X.SampleID" in d.columns else d.columns[0]
    d["sample"] = d[key].astype(str).str.extract(_SAMPLE, expand=False)
    d = d.dropna(subset=["sample"])
    num = d.select_dtypes("number").columns
    agg = {c: "mean" for c in num}
    agg["UserName"] = "first"
    if "StudyDayNo" in d.columns:
        agg["StudyDayNo"] = "first"
    g = d.groupby("sample").agg(agg)

    kcal = g["KCAL"].replace(0, np.nan)
    out = pd.DataFrame(index=g.index)
    out["person"] = g["UserName"].astype(str)
    out["day"] = (g["StudyDayNo"].astype(str).str.extract(r"(\d+)", expand=False).astype(float)
                  if "StudyDayNo" in g.columns else np.nan)
    out["kcal"] = g["KCAL"]
    out["fiber"] = g.get("FIBE")
    out["fiber_dens"] = g.get("FIBE") / kcal * 1000.0           # g per 1000 kcal (energy-adjusted)

    def _z(col):
        s = (g[col] / kcal) if col in g.columns else pd.Series(np.nan, index=g.index)
        return (s - s.mean()) / s.std(ddof=0) if s.std(ddof=0) else s * 0.0
    good = [c for c in _DQ_GOOD if c in g.columns]
    bad = [c for c in _DQ_BAD if c in g.columns]
    out["dq_score"] = pd.concat([_z(c) for c in good], axis=1).mean(axis=1) \
        - pd.concat([_z(c) for c in bad], axis=1).mean(axis=1)
    return out


def build(repo: Path, out_path: Path) -> pd.DataFrame:
    enzyme = repo / "data" / "function" / "KEGG" / "kegg_enzyme_table.txt"
    diet = repo / "data" / "diet" / "nutrition_totals.txt"
    for p in (enzyme, diet):
        if not p.exists():
            raise SystemExit(f"missing {p}\nClone the public repo first:\n  git clone "
                             f"https://github.com/knights-lab/dietstudy_analyses {repo}")
    axes = _ec_axes(enzyme)
    dv = _diet_vars(diet)
    df = dv.join(axes, how="inner")                            # shared MCT.f.#### key
    df = df.reset_index().rename(columns={"index": "sample"})
    cols = ["sample", "person", "day", "kcal", "fiber", "fiber_dens", "dq_score",
            "estrobolome", "glycogen_starch", "scfa"]
    df = df[[c for c in cols if c in df.columns]]
    df.to_csv(out_path, sep="\t", index=False)
    print(f"wrote {out_path}  ({df['person'].nunique()} people, {len(df)} person-days, "
          f"axes: {[c for c in ('estrobolome','glycogen_starch','scfa') if c in df.columns]})")
    return df


def main(argv=None):
    p = argparse.ArgumentParser(description="build within-person diet->function input from Johnson 2019")
    p.add_argument("--repo", default=str(HERE / "data" / "dietstudy_analyses"),
                   help="local clone of knights-lab/dietstudy_analyses")
    p.add_argument("--out", default=str(HERE / "data" / "diet_function_input.tsv"))
    a = p.parse_args(argv)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    build(Path(a.repo), Path(a.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
