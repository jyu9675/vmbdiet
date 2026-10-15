"""Generate an Aim-3 functional demo: a HUMAnN 3 EC-level merged abundance table across a
repeated-measures gut-vagina cohort, with planted diet->function signals, plus matching metadata.

EC-level (not pathway-level) so all three axis categories are present in one table:
  * estrobolome      3.2.1.31  beta-glucuronidase          -> rises with ANIMAL_PROTEIN
  * glycogen_starch  3.2.1.1   alpha-amylase               -> rises with FIBER
                     3.2.1.3   glucan 1,4-alpha-glucosidase -> rises with FIBER
  * scfa             2.7.2.7   butyrate kinase             -> rises with FIBER
                     2.7.2.1   acetate kinase              -> rises with FIBER
                     2.8.3.8   butyryl-CoA:acetate CoA-transferase -> flat (null)
Plus null ECs and the HUMAnN special rows. Regroup real data with
`humann_regroup_table --groups uniref90_level4ec` to get this layout.
"""
import numpy as np
import pandas as pd
from pathlib import Path

EC_ROWS = {
    "3.2.1.31: beta-glucuronidase": ("animal_protein", 0.012),
    "3.2.1.1: alpha-amylase": ("fiber", 0.05),
    "3.2.1.3: glucan 1,4-alpha-glucosidase": ("fiber", 0.04),
    "2.7.2.7: butyrate kinase": ("fiber", 0.045),
    "2.7.2.1: acetate kinase": ("fiber", 0.03),
    "2.8.3.8: acetate CoA-transferase": (None, 0.0),       # SCFA gene, no diet effect (null control)
    "1.1.1.1: alcohol dehydrogenase": (None, 0.0),         # null filler
    "2.7.7.7: DNA-directed DNA polymerase": (None, 0.0),   # null filler
    "6.3.5.5: carbamoyl-phosphate synthase": (None, 0.0),  # null filler
}


def main(n_part=48, timepoints=2, seed=13):
    rng = np.random.default_rng(seed)
    fiber = rng.normal(25, 8, n_part)
    animal_protein = rng.normal(60, 18, n_part)
    samples, meta = [], []
    for p in range(n_part):
        for site in ("rectal", "vaginal"):
            for t in range(timepoints):
                samples.append(f"P{p:03d}_{site[:3]}_t{t}")
                meta.append(dict(participant=f"P{p:03d}", body_site=site, timepoint=t,
                                 fiber=round(fiber[p] + rng.normal(0, 1), 2),
                                 animal_protein=round(animal_protein[p] + rng.normal(0, 2), 2),
                                 age=round(float(rng.normal(30, 5)), 1)))
    meta = pd.DataFrame(meta, index=samples)

    rows = {}
    for ec, (driver, beta) in EC_ROWS.items():
        vals = []
        for s in samples:
            base = 2.0
            if driver:
                x = meta.loc[s, driver]
                center = 25 if driver == "fiber" else 60
                base += beta * (x - center)
            vals.append(np.exp(base + rng.normal(0, 0.3)) * 1e4)   # CPM-ish scale
        rows[ec] = vals
    feat = pd.DataFrame(rows, index=samples).T                      # features x samples (HUMAnN layout)
    # HUMAnN special rows
    feat.loc["UNMAPPED"] = rng.uniform(1e5, 2e5, len(samples))
    feat.loc["UNGROUPED"] = rng.uniform(2e5, 4e5, len(samples))
    feat.index.name = "# Gene Family"

    out = Path(__file__).parent / "data"
    out.mkdir(exist_ok=True)
    feat.to_csv(out / "ecs_merged.tsv", sep="\t")
    meta.to_csv(out / "metadata.tsv", sep="\t")
    print(f"wrote {out/'ecs_merged.tsv'} ({feat.shape[0]} features x {feat.shape[1]} samples) and metadata.tsv")


if __name__ == "__main__":
    main()
