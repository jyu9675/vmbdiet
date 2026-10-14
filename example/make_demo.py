#!/usr/bin/env python3
"""Generate a synthetic demo (abundance + diet metadata) with a built-in diet->microbiome signal,
so `vmbdiet associate` has something to run on. NOT real data — illustrative only."""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
TAXA = ["Lactobacillus_crispatus", "Lactobacillus_iners", "Lactobacillus_gasseri",
        "Gardnerella_vaginalis", "Prevotella_bivia", "Atopobium_vaginae", "Sneathia_amnii"]

def main(n=200, seed=7):
    rng = np.random.default_rng(seed)
    carb = rng.normal(250, 60, n); fat = rng.normal(70, 20, n)
    folate = rng.normal(400, 120, n); age = rng.normal(30, 5, n)
    race = rng.choice(["White", "Black", "Hispanic", "Asian"], n)
    bmi = rng.normal(26, 5, n)
    # crispatus propensity: up with carb & folate, down with fat
    z = 0.018*(carb-250) + 0.004*(folate-400) - 0.03*(fat-70) + rng.normal(0, 0.5, n)
    rows = []
    for i in range(n):
        if z[i] > 0.3:
            comp = [rng.uniform(0.7, 0.95), 0.05, 0.03, 0.02, 0.02, 0.02, 0.01]
        elif z[i] < -0.3:
            comp = [0.03, 0.05, 0.02, rng.uniform(0.3, 0.5), rng.uniform(0.15, 0.3), 0.1, 0.08]
        else:
            comp = [0.05, rng.uniform(0.6, 0.85), 0.04, 0.05, 0.03, 0.02, 0.02]
        comp = np.array(comp); comp = comp / comp.sum()
        rows.append((comp * rng.integers(5000, 25000)).round().astype(int))
    idx = [f"V{i:03d}" for i in range(n)]
    pd.DataFrame(rows, columns=TAXA, index=idx).to_csv(os.path.join(HERE, "demo_abund.tsv"), sep="\t")
    pd.DataFrame({"carbohydrate_g": carb.round(1), "fat_g": fat.round(1),
                  "folate_ug": folate.round(1), "age": age.round(1),
                  "bmi": bmi.round(1), "race": race}, index=idx).to_csv(
                  os.path.join(HERE, "demo_meta.tsv"), sep="\t")
    print(f"wrote demo_abund.tsv ({n}×{len(TAXA)}) + demo_meta.tsv")

if __name__ == "__main__":
    main()
