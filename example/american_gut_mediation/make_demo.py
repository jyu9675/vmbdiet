"""Synthetic stand-in for the American Gut mediation input, to validate the harness offline (no
download). Shape mirrors fetch_americangut.py's output: one row per participant with a diet-quality
score (HEI), gut alpha diversity (shannon), BMI, age, sex.

Planted causal structure (a PARTIAL mediation): better diet -> higher gut diversity -> lower BMI,
plus a direct diet->BMI path. A correct harness recovers a negative, significant indirect effect
(NIE) through shannon and a proportion mediated strictly between 0 and 1.
"""
import numpy as np
import pandas as pd
from pathlib import Path


def make(n=600, seed=7):
    rng = np.random.default_rng(seed)
    hei = rng.normal(60, 12, n).clip(20, 100)                 # HEI-2010 diet quality
    age = rng.normal(45, 14, n).clip(18, 90)
    sex = rng.choice(["female", "male"], n)
    # mediator: diversity rises with diet quality        (path a > 0)
    shannon = 3.5 + 0.015 * (hei - 60) + rng.normal(0, 0.5, n)
    # outcome: BMI falls with diversity (path b < 0) AND directly with diet (c' < 0)
    bmi = (27.0 - 0.030 * (hei - 60) - 1.2 * (shannon - 3.5)
           + rng.normal(0, 3.0, n)).clip(14, 55)
    df = pd.DataFrame(dict(sample_id=[f"demo.{i:04d}" for i in range(n)],
                           hei=hei, shannon=shannon, bmi=bmi, age=age, sex=sex,
                           healthy_subset="No"))
    return df


def main():
    out = Path(__file__).parent / "data"
    out.mkdir(exist_ok=True)
    df = make()
    path = out / "agp_demo_input.tsv"
    df.to_csv(path, sep="\t", index=False)
    print(f"wrote {path} ({len(df)} synthetic participants)")
    return str(path)


if __name__ == "__main__":
    main()
