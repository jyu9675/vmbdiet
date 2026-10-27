#!/usr/bin/env python3
"""Fetch the public American Gut manuscript-package mapping file and distill it to a compact,
analysis-ready TSV for the diet -> gut -> health mediation prototype. Stdlib only (urllib).

American Gut / The Microsetta Initiative released, per rarefaction depth, a single mapping file
`ag_map_with_alpha.txt` that already carries (a) the VioScreen FFQ results appended per participant,
(b) alpha diversity (Shannon, Faith's PD, observed OTUs) at that depth and lower, and (c) the AGP
survey metadata (BMI, age, sex, and a `healthy_subset` flag). One fecal sample per participant. So a
fully-public diet+microbiome+health table needs just this one download — no BIOM crunching, no login.

Source (CC-licensed public data): https://ftp.microbio.me/AmericanGut/manuscript-package/<depth>/
We keep only the columns the mediation needs and write them to `agp_mediation_input.tsv`:
  sample_id, hei (diet quality), shannon (gut alpha diversity), bmi, age, sex, healthy_subset

Usage:
  python fetch_americangut.py                 # depth 10000 (manuscript primary), -> agp_mediation_input.tsv
  python fetch_americangut.py --depth 2500    # shallower rarefaction keeps more samples
  python fetch_americangut.py --map local_ag_map.txt   # use an already-downloaded map file
"""
import argparse
import csv
import io
import os
import sys
import urllib.request
from pathlib import Path

HERE = Path(__file__).parent
BASE = "https://ftp.microbio.me/AmericanGut/manuscript-package"
MISSING = {"", "NA", "no_data", "Not provided", "Unspecified", "Unknown", "not applicable",
           "None", "nan", "Missing: Not provided"}


def _num(v):
    if v is None or str(v).strip() in MISSING:
        return ""
    try:
        return float(str(v).strip())
    except ValueError:
        return ""


def main(argv=None):
    p = argparse.ArgumentParser(description="distill the AGP manuscript map to a mediation input TSV")
    p.add_argument("--depth", default="10000", help="rarefaction depth folder (1250/2500/5000/10000)")
    p.add_argument("--map", help="path to an already-downloaded ag_map_with_alpha.txt (skip download)")
    p.add_argument("--out", default=str(HERE / "agp_mediation_input.tsv"))
    a = p.parse_args(argv)

    hei_c, bmi_c, age_c, sex_c, flag_c = ("vioscreen_hei2010_score", "bmi_corrected",
                                          "age_corrected", "sex", "healthy_subset")
    sh_c = f"shannon_{a.depth}"

    if a.map:
        txt = Path(a.map).read_text(encoding="utf-8", errors="replace")
    else:
        url = f"{BASE}/{a.depth}/ag_map_with_alpha.txt"
        print(f"downloading {url} ...")
        txt = urllib.request.urlopen(url, timeout=180).read().decode("utf-8", "replace")

    rd = csv.DictReader(io.StringIO(txt), delimiter="\t")
    idcol = rd.fieldnames[0]                                   # '#SampleID'
    for c in (hei_c, sh_c, bmi_c, age_c, sex_c):
        if c not in rd.fieldnames:
            sys.exit(f"column {c!r} not in map (depth {a.depth}?) — available shannon cols: "
                     + ", ".join(c for c in rd.fieldnames if c.startswith('shannon_')))

    kept = 0
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "hei", "shannon", "bmi", "age", "sex", "healthy_subset"])
        for r in rd:
            hei, sh, bmi, age = _num(r[hei_c]), _num(r[sh_c]), _num(r[bmi_c]), _num(r[age_c])
            sex = (r.get(sex_c) or "").strip().lower()
            if "" in (hei, sh, bmi, age) or sex not in ("male", "female"):
                continue                                       # need the full chain + clean sex
            w.writerow([r[idcol], hei, sh, bmi, age, sex, (r.get(flag_c) or "").strip()])
            kept += 1
    print(f"wrote {a.out}: {kept} participants with HEI + shannon_{a.depth} + BMI + age + sex")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
