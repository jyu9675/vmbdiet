"""Generate a gut-vagina-axis mediation demo: does dietary fiber lower CST-IV risk THROUGH gut SCFA
capacity? Plants a partial mediation so `vmbdiet mediate` has a known truth to recover.

Causal model built in (one sample per participant — the clean cross-sectional case):
    fiber (X)  -->  gut_scfa capacity (M)  -->  CST-IV (Y, binary)
    - fiber raises gut SCFA capacity            (a > 0)
    - higher gut SCFA lowers CST-IV risk        (b < 0, the mediated path — the bulk of the effect)
    - fiber also has a small DIRECT protective effect not through SCFA  (c' < 0, small)
So most of fiber's protection against CST-IV should flow through the gut (proportion mediated ~0.6-0.8).

In a real run the mediator column is a `vmbdiet functional` axis capacity score (gut-site) and the
outcome is a `vmbdiet cst` call on the vaginal sample; here both are simulated.
"""
import numpy as np
import pandas as pd
from pathlib import Path


def main(n=300, seed=17):
    rng = np.random.default_rng(seed)
    fiber = rng.normal(25, 8, n)
    age = rng.normal(30, 5, n)
    gut_scfa = 2.0 + 0.05 * (fiber - 25) + rng.normal(0, 0.3, n)          # mediator (log-capacity-ish)
    zf = (fiber - fiber.mean()) / fiber.std()
    zs = (gut_scfa - gut_scfa.mean()) / gut_scfa.std()
    za = (age - age.mean()) / age.std()
    # CST-IV risk: strongly lowered by gut SCFA (mediated), weakly by fiber directly
    logit = -0.4 - 1.3 * zs - 0.3 * zf + 0.2 * za
    cst_IV = (rng.random(n) < 1 / (1 + np.exp(-logit))).astype(int)
    df = pd.DataFrame({"participant": [f"P{i:03d}" for i in range(n)],
                       "fiber": fiber.round(2), "age": age.round(1),
                       "gut_scfa": gut_scfa.round(3), "cst_IV": cst_IV},
                      index=[f"S{i:03d}" for i in range(n)])
    out = Path(__file__).parent / "data"; out.mkdir(exist_ok=True)
    df.to_csv(out / "axis_table.tsv", sep="\t")
    print(f"wrote {out/'axis_table.tsv'}  (n={n}, CST-IV prevalence={cst_IV.mean():.2f})")


if __name__ == "__main__":
    main()
