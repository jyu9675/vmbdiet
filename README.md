# vmbdiet

**Diet ↔ vaginal microbiome association toolkit.** Does what a person eats track with their vaginal community — bacterial vaginosis (BV), community state type (CST), *Lactobacillus* dominance? This packages the statistical methods shared across the diet–vaginal-microbiome literature into one tested, reusable tool, so the association can be run on any cohort that pairs a 16S table with dietary data.

> Status **v0.1** — engine, CLI, and synthetic validation complete (recovers injected diet→microbiome signals). Real public-dataset validation in progress; see [docs/DATASETS.md](docs/DATASETS.md).

## Why

Four studies converge on a real signal, each with different data and methods:

| Study | Finding |
|---|---|
| Neggers 2007 (BV) | ↑ dietary **fat** → ↑ BV risk; ↑ **folate / vitamin A / calcium** → ↓ BV risk |
| Tuddenham 2019 (prenatal) | ↑ **dairy / fruit / fiber / vitamin D** → *L. crispatus* over *L. iners* (dose-dependent) |
| Hawaii multi-ethnic 2024 | ↑ **carbohydrate** → *L. crispatus*; low carb → *L. iners* |
| Italy 2025 (public data) | ↑ **animal protein** → CST IV; **alcohol** → *Gardnerella* |

They all run the same core analysis — classify the community, then regress/correlate diet against it, adjusting for covariates. `vmbdiet` is that core, standardized.

## What it does

- **CST classification** ([cst.py](vmbdiet/cst.py)) — assign each 16S sample a community state type (I = *L. crispatus*, II = *L. gasseri*, III = *L. iners*, V = *L. jensenii*, IV = diverse/BV-associated). A transparent dominant-taxon classifier by default; VALENCIA centroid + Yue-Clayton scoring when you supply the reference centroids.
- **Adjusted association** ([associate.py](vmbdiet/associate.py)) — logistic regression of an outcome (BV / CST-IV / *Lactobacillus*-depleted) on each nutrient **adjusted for covariates** (age, race, BMI, …), as per-SD odds ratios with 95% CI and BH-FDR; nutrient↔taxon correlations; Mann-Whitney between dominant/depleted groups; **PERMANOVA** of community composition against a grouping.
- **Diversity** ([diversity.py](vmbdiet/diversity.py)) — Bray-Curtis, Shannon, PCoA for the beta-diversity + PERMANOVA workflow.

## Install

```bash
pip install -e .          # deps: pandas, numpy, scipy, statsmodels
```

## Quickstart (shipped demo)

```bash
python example/make_demo.py          # writes a synthetic abundance + diet table
vmbdiet associate --abund example/demo_abund.tsv --meta example/demo_meta.tsv \
    --outcome cst_IV --nutrients carbohydrate_g,fat_g,folate_ug --covariates age,race,bmi \
    --out example/demo_out
```

```
PERMANOVA community~cst_IV: pseudo-F=320.24, p=0.001

adjusted logistic (cst_IV ~ nutrient + 3 covariates), per-SD OR:
  nutrient                   OR           95% CI        p        q
  carbohydrate_g           0.14        0.08-0.25    0.000    0.000 *
  fat_g                    2.17        1.54-3.06    0.000    0.000 *
  folate_ug                0.54        0.39-0.75    0.000    0.000 *
```

The synthetic data has fat injected as a risk factor and carbohydrate/folate as protective — and the tool recovers exactly that (fat OR > 1, carb/folate OR < 1), the same directions the literature reports. (Illustrative data, not a real cohort.)

## Inputs

1. **Abundance table** — TSV, samples (rows) × taxa (columns), counts or relative abundance. Taxa named so species are recognizable (`Lactobacillus_crispatus`, `Gardnerella_vaginalis`, …).
2. **Metadata table** — TSV, samples × (nutrient/diet columns + covariates [+ an outcome column]). `vmbdiet` also derives `cst_IV` and `lacto_depleted` from the abundance, so you can model those directly.

## Outputs

`cst.tsv` (per-sample CST), `logistic_scan.tsv` (adjusted OR/CI/p/q per nutrient), `nutrient_taxon_corr.tsv`, `summary.json` (CST distribution, PERMANOVA, top hits).

## Relationship to the rest of the program

Sibling to the strain-level tools (strainshare, engrafttracker). Those ask *which strain, and did it transfer/engraft*; `vmbdiet` asks *does host diet shape which community establishes* — the modifiable-exposure side of the same vaginal-health question, and the regime (low-diversity, 16S) that's computationally tractable on a laptop.

## License

MIT © Jean Yu
