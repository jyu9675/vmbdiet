"""Synthetic stand-in for several curatedMetagenomicData cohorts, to validate the cross-cohort
function <-> cardiometabolic reproducibility harness offline (no R, no cMD download).

Each cohort is a HUMAnN 3 pathway table (features x samples) + a cMD-style metadata table carrying a
cardiometabolic `disease` label (Healthy / Type 2 Diabetes Mellitus), age, sex, BMI, westernized=Yes
(cMD's cardiometabolic cohorts are uniformly westernized -- see docs/diet-gut-health-landscape.md).
The planted signal is reproducible across all three cohorts: **SCFA (butyrate) biosynthesis capacity
is lower in T2D** (the well-known butyrate-depletion direction), glycogen/starch is null. A correct
harness recovers a consistent pooled OR < 1 for SCFA (low I^2) and a non-significant glycogen pool.
Estrobolome (beta-glucuronidase) is gene-level (EC/KO), not a pathway, so it does not appear here.
"""
import numpy as np
import pandas as pd
from pathlib import Path

# MetaCyc pathway IDs the axis dictionary recognizes (+ nulls); (axis, log-shift applied to T2D cases)
PATHWAYS = {
    "GLYCOCAT-PWY: glycogen degradation I": ("glycogen", 0.0),
    "PWY-6737: starch degradation V": ("glycogen", 0.0),
    "PWY-5676: acetyl-CoA fermentation to butanoate II": ("scfa", -0.30),   # modestly lower SCFA in T2D
    "CENTFERM-PWY: pyruvate fermentation to butanoate": ("scfa", -0.25),
    "PWY-7219: adenosine ribonucleotides de novo biosynthesis": (None, 0.0),
    "PWY-6700: queuosine biosynthesis": (None, 0.0),
}

# (study_name, n_healthy, n_t2d, seed) -- three independent cohorts of different sizes, like cMD's
# KarlssonFH_2013 / QinJ_2012 / SankaranarayananK_2015.
COHORTS = [("CohortA_2013", 60, 55, 11), ("CohortB_2012", 80, 70, 23), ("CohortC_2015", 24, 24, 37)]


def one(study, n_ctrl, n_case, seed):
    rng = np.random.default_rng(seed)
    samples, meta = [], []
    for i in range(n_ctrl + n_case):
        case = 1 if i >= n_ctrl else 0
        sid = f"{study}_S{i:03d}"
        samples.append(sid)
        meta.append(dict(sample_id=sid, study_name=study, T2D=case,
                         disease=("Type 2 Diabetes Mellitus" if case else "Healthy"),
                         age=round(float(rng.normal(60 if case else 58, 8)), 1),
                         sex=str(rng.choice(["male", "female"])),
                         BMI=round(float(rng.normal(29 if case else 26, 4)), 1),
                         westernized="Yes"))
    md = pd.DataFrame(meta).set_index("sample_id")
    rows = {}
    for pw, (driver, beta) in PATHWAYS.items():
        vals = []
        for sid in samples:
            case = int(md.loc[sid, "T2D"])
            base = 2.0 + (beta if (driver and case) else 0.0)
            vals.append(np.exp(base + rng.normal(0, 0.6)) * 1e3)   # CPM-ish (wide within-group spread)
        rows[pw] = vals
    feat = pd.DataFrame(rows, index=samples).T
    feat.loc["UNMAPPED"] = rng.uniform(1e4, 2e4, len(samples))
    feat.loc["UNINTEGRATED"] = rng.uniform(2e4, 4e4, len(samples))
    feat.index.name = "# Pathway"
    return feat, md


def main():
    out = Path(__file__).parent / "data"
    out.mkdir(exist_ok=True)
    cohorts = []
    for study, nc, na, seed in COHORTS:
        feat, md = one(study, nc, na, seed)
        p = out / f"{study}_pathway_abundance.tsv"
        m = out / f"{study}_metadata.tsv"
        feat.to_csv(p, sep="\t")
        md.to_csv(m, sep="\t")
        cohorts.append((study, str(p), str(m)))
        print(f"wrote {study}: {feat.shape[0]}x{feat.shape[1]} pathways, {nc} Healthy / {na} T2D")
    return cohorts


if __name__ == "__main__":
    main()
