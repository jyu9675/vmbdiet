"""Generate a small gut-vagina-axis demo for the MaAsLin 3 bridge and write MaAsLin-format TSVs.

Shape mirrors the FRESH design: each participant contributes rectal + vaginal samples across two
timepoints. Diet ('fiber', 'animal_protein') is per-participant (FFQ-style). Three features carry a
planted biology so you can see what `vmbdiet maaslin` recovers:
  * PWY_glycogen_degradation  — abundance RISES with fiber     (estrobolome/glycogen mechanism stand-in)
  * GUS_beta_glucuronidase    — PREVALENCE rises with fiber    (estrobolome gene family stand-in)
  * SCFA_butyrate_synth       — abundance RISES with animal_protein
The rest are null/filler (one driven only by participant identity, to show the random effect at work).
"""
import numpy as np
import pandas as pd
from pathlib import Path

FEATURES = ["PWY_glycogen_degradation", "GUS_beta_glucuronidase", "SCFA_butyrate_synth",
            "PWY_participant_null", "PWY_flat", "PWY_noise"]


def main(n_part=50, timepoints=2, seed=11):
    rng = np.random.default_rng(seed)
    fiber = rng.normal(25, 8, n_part)
    animal_protein = rng.normal(60, 18, n_part)
    base = rng.normal(0, 1.3, n_part)          # participant intercept for the null feature
    rows, meta, index = [], [], []
    for p in range(n_part):
        for site in ("rectal", "vaginal"):
            for t in range(timepoints):
                fb = fiber[p] + rng.normal(0, 1.0)
                ap = animal_protein[p] + rng.normal(0, 2.0)
                glyc = np.exp(2.0 + 0.07 * (fb - 25) + rng.normal(0, 0.3))
                gus_p = 1 / (1 + np.exp(-(0.22 * (fb - 25))))
                gus = (rng.random() < gus_p) * np.exp(1.4 + rng.normal(0, 0.3))
                scfa = np.exp(2.0 + 0.03 * (ap - 60) + rng.normal(0, 0.3))
                null = np.exp(2.0 + 1.4 * base[p] + rng.normal(0, 0.25))
                flat = np.exp(3.0 + rng.normal(0, 0.2))
                noise = np.exp(1.5 + rng.normal(0, 0.8))
                rows.append([glyc, gus, scfa, null, flat, noise])
                meta.append(dict(participant=f"P{p:03d}", body_site=site, timepoint=t,
                                 fiber=round(fb, 2), animal_protein=round(ap, 2),
                                 age=round(float(rng.normal(30, 5)), 1)))
                index.append(f"P{p:03d}_{site[:3]}_t{t}")
    out = Path(__file__).parent / "data"
    out.mkdir(exist_ok=True)
    pd.DataFrame(rows, columns=FEATURES, index=index).to_csv(out / "features.tsv", sep="\t")
    pd.DataFrame(meta, index=index).to_csv(out / "metadata.tsv", sep="\t")
    print(f"wrote {out/'features.tsv'} and {out/'metadata.tsv'}  ({len(index)} samples, {n_part} participants)")


if __name__ == "__main__":
    main()
