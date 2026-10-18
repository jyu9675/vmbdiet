"""Synthetic stand-in for a curatedMetagenomicData pathway export, to validate the westernized-contrast
Aim-3 harness offline (no R, no cMD download).

Shape mirrors what `write_cmd_export_R(..., data_type="pathway_abundance")` produces: a HUMAnN 3 pathway
table (features x samples) + a cMD-style metadata table carrying `westernized`. Westernization is cMD's
strongest diet/lifestyle axis (non-westernized ~ high-fiber, plant-rich), so the planted signal is
non-westernized -> higher glycogen-degradation and SCFA-biosynthesis capacity (the fiber-fermentation
expectation). Estrobolome (beta-glucuronidase) is gene-level (EC/KO), not a MetaCyc pathway, so it needs
the gene-family export instead and does not appear here.
"""
import numpy as np
import pandas as pd
from pathlib import Path

# MetaCyc pathway IDs the axis dictionary recognizes (+ nulls); (axis driver, effect in non-westernized)
PATHWAYS = {
    "GLYCOCAT-PWY: glycogen degradation I": ("glycogen", +0.6),
    "PWY-6737: starch degradation V": ("glycogen", +0.5),
    "PWY-5676: acetyl-CoA fermentation to butanoate II": ("scfa", +0.6),
    "CENTFERM-PWY: pyruvate fermentation to butanoate": ("scfa", +0.5),
    "PWY-7219: adenosine ribonucleotides de novo biosynthesis": (None, 0.0),
    "PWY-6700: queuosine biosynthesis": (None, 0.0),
    "DTDPRHAMSYN-PWY: dTDP-L-rhamnose biosynthesis I": (None, 0.0),
}


def main(n_west=40, n_nonwest=25, seed=19):
    rng = np.random.default_rng(seed)
    samples, meta = [], []
    for i in range(n_west + n_nonwest):
        west = "Westernized" if i < n_west else "Non-westernized"
        samples.append(f"S{i:03d}")
        meta.append(dict(sample_id=f"S{i:03d}", study_name="DemoWesternized_2026",
                         westernized=west, country=("USA" if west == "Westernized" else "FJI"),
                         age=round(float(rng.normal(35, 8)), 1)))
    meta = pd.DataFrame(meta).set_index("sample_id")

    rows = {}
    for pw, (driver, beta) in PATHWAYS.items():
        vals = []
        for s in samples:
            base = 2.0
            if driver and meta.loc[s, "westernized"] == "Non-westernized":
                base += beta
            vals.append(np.exp(base + rng.normal(0, 0.3)) * 1e3)   # CPM-ish
        rows[pw] = vals
    feat = pd.DataFrame(rows, index=samples).T
    feat.loc["UNMAPPED"] = rng.uniform(1e4, 2e4, len(samples))
    feat.loc["UNINTEGRATED"] = rng.uniform(2e4, 4e4, len(samples))
    feat.index.name = "# Pathway"

    out = Path(__file__).parent / "data"
    out.mkdir(exist_ok=True)
    feat.to_csv(out / "cmd_demo_pathway_abundance.tsv", sep="\t")
    meta.to_csv(out / "cmd_demo_metadata.tsv", sep="\t")
    print(f"wrote {out/'cmd_demo_pathway_abundance.tsv'} ({feat.shape[0]}x{feat.shape[1]}) and metadata "
          f"({n_west} Westernized / {n_nonwest} Non-westernized)")


if __name__ == "__main__":
    main()
