#!/usr/bin/env python3
"""Synthetic stand-in for the Johnson 2019 within-person diet<->function data, so analyze.py and the
test run with no download. Plants a PERSONALIZED (heterogeneous) within-person fiber -> SCFA-capacity
coupling that is positive in most people (pooled r>0, significant, 'consistent' YES, with real I^2 from
between-person slope differences), a NULL glycogen axis (slopes centered on 0 -> not consistent), and a
modest estrobolome signal. Numbers are synthetic — a pipeline check, not a finding."""
import numpy as np
import pandas as pd

N_PEOPLE = 24
N_DAYS = 15


def make(seed: int = 0) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    # per-person within-person slopes (standardized units): scfa mostly +, glycogen null, estrobolome +.
    scfa_b = rng.normal(0.60, 0.45, N_PEOPLE)             # heterogeneous but mostly positive
    glyc_b = rng.normal(0.00, 0.45, N_PEOPLE)             # null: ~half each sign
    estr_b = rng.normal(0.45, 0.55, N_PEOPLE)             # modest positive
    rows = []
    for p in range(N_PEOPLE):
        person = f"S{p+1:02d}"
        kcal = rng.normal(2100, 350, N_DAYS).clip(1200, 3600)
        fiber_dens = rng.normal(12.0, 4.0, N_DAYS).clip(3, 30)      # g / 1000 kcal
        fz = (fiber_dens - fiber_dens.mean()) / (fiber_dens.std() or 1.0)
        dq = fz * 0.7 + rng.normal(0, 0.7, N_DAYS)                  # diet-quality tracks fiber loosely
        scfa = 0.020 + 0.004 * scfa_b[p] * fz + rng.normal(0, 0.004, N_DAYS)
        glyc = 0.015 + 0.003 * glyc_b[p] * fz + rng.normal(0, 0.004, N_DAYS)
        estr = 0.0020 + 0.0005 * estr_b[p] * fz + rng.normal(0, 0.0006, N_DAYS)
        for d in range(N_DAYS):
            rows.append(dict(sample=f"{person}.d{d+1:02d}", person=person, day=d + 1,
                             kcal=kcal[d], fiber=fiber_dens[d] * kcal[d] / 1000.0,
                             fiber_dens=fiber_dens[d], dq_score=dq[d],
                             estrobolome=max(estr[d], 1e-6), glycogen_starch=max(glyc[d], 1e-6),
                             scfa=max(scfa[d], 1e-6)))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    df = make()
    print(df.head())
    print(f"\n{df['person'].nunique()} people x ~{N_DAYS} days = {len(df)} person-days")
