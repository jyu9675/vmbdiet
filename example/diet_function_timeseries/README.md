# Within-person daily diet → gut FUNCTION (Johnson 2019 "Daily Sampling")

The [American Gut mediation](../american_gut_mediation/) left one question open: diet→BMI was real, but
gut **α-diversity** mediated ~none of it — was that because diversity is too coarse a mediator, or
because gut function genuinely doesn't carry the signal? Answering it needs *measured functional
capacity* paired with diet, which the landscape doc said public data lacked. This example uses the one
fully-public dataset that has it: **Johnson et al. 2019, *Cell Host & Microbe*** — 34 healthy adults,
**17 consecutive days** of daily **ASA24 food records + daily shotgun metagenomes**
([knights-lab/dietstudy_analyses](https://github.com/knights-lab/dietstudy_analyses)).

It is a **within-person** design, not a BMI-outcome cohort: BMI and blood markers are measured once per
person (n≈34), so this tests **diet → function dynamics**, not diet → BMI. What it uniquely enables is
the **estrobolome gene-level read** — β-glucuronidase **EC 3.2.1.31** — which curatedMetagenomicData
could not deliver cheaply (its `gene_families` export is ~6M ID-only UniRef90 rows). Here the shotgun
functional tables are already KEGG-annotated, so the existing `vmbdiet.functional` **EC** dictionaries
match all three axes directly (no MetaCyc-pathway regrouping needed).

## Method — the "unit" is a PERSON

Mirrors the [cross-cohort harness](../gut_function_cardiometabolic/), but each **person** is their own
control:

1. for each person, correlate a daily diet variable with a daily axis capacity across their ~17 days
   (within-person Pearson *r*) — stable between-person differences in body size, transit, enterotype
   **cannot** confound a within-person correlation;
2. pool the per-person correlations across people by **Fisher-z inverse-variance meta-analysis**
   (`vmbdiet.meta_analyze`, DerSimonian–Laird) → a pooled within-person effect, CI, and **I²**. High I²
   would itself be the original paper's finding: diet–microbiome coupling is *personalized*.

Axes (`vmbdiet.functional`, matched on KEGG enzyme EC numbers): **estrobolome** (β-glucuronidase, EC
3.2.1.31), **scfa** (fermentation enzymes — the fiber mechanism), **glycogen_starch** (amylases).

## Run it

Real data is **not redistributed** (gitignored); reproduce from the public repo:

```bash
git clone https://github.com/knights-lab/dietstudy_analyses \
  example/diet_function_timeseries/data/dietstudy_analyses
python fetch_dietstudy.py                                   # -> data/diet_function_input.tsv
python analyze.py --data data/diet_function_input.tsv --diet-var fiber_dens   # also: fiber, dq_score
python analyze.py --data data/diet_function_input.tsv --diet-var fiber_dens --lag 1   # prior-day diet
```

Offline (no download) runs the synthetic demo, which plants a personalized fiber→SCFA coupling the
harness recovers (locked by `tests/test_diet_function_timeseries.py`):

```bash
python analyze.py          # SYNTHETIC: scfa consistent YES (pooled r≈+0.51), glycogen no
```

## Real-data result (37 people, 521 paired person-days) — a reproducible null

Diet variable = dietary-fiber density (g / 1000 kcal), same-day pairing:

| axis | k people | pooled within-person r [95% CI] | p | I² | % same-sign | consistent |
|---|---|---|---|---|---|---|
| estrobolome (β-glucuronidase) | 33 | +0.086 [−0.03, +0.20] | 0.135 | 11% | 61% | **no** |
| scfa | 33 | −0.042 [−0.15, +0.06] | 0.439 | 0% | 61% | **no** |
| glycogen_starch | 33 | +0.003 [−0.10, +0.11] | 0.956 | 0% | 52% | **no** |

The null is **stable across diet variables (absolute fiber, fiber density, a plant/fiber-forward
diet-quality score) and across a 1-day diet→stool lag.** One borderline hit appeared — diet-quality
score → glycogen capacity, r −0.106, p=0.050 — but it sits exactly at the threshold, is one of ~24
configurations tested, and points opposite to the mechanistic prior, so it does **not** survive
multiple comparisons and is reported as noise, not a finding.

**Positive controls confirm the null is real, not a mis-join:** the same within-person machinery finds
fiber ↔ energy intake at median r +0.58 and fiber ↔ fiber-density at +0.72, and the axis capacities do
vary day to day (within-person CV 0.13–0.26). The pipeline detects coupling where it exists.

## Conclusion

**Daily diet does not reproducibly move coarse, summed EC-level functional capacity** (SCFA,
estrobolome, glycogen) within person in this cohort — and with **I² ≈ 0**, it is not even a
personalization story at this resolution; the summed-capacity readout simply doesn't track daily diet.
This is coherent with the original paper's deeper point: diet–microbiome associations are personalized
**and** live at *fine* feature resolution (a specific food ↔ a specific microbe/gene), so a broad
"total fiber → total SCFA-enzyme capacity" axis washes out.

For the gut–vagina program this is a useful negative: it says the estrobolome / SCFA **axis-capacity**
abstraction — the summed-pathway proxy used throughout `vmbdiet.functional` — is too coarse to detect a
diet response even with *daily shotgun* data and a person as their own control. A real test of the
estrobolome mechanism needs finer resolution: the GUS **structural subclass** that actually reactivates
estrogens (Ervin 2019), not total β-glucuronidase, and ideally measured SCFA/estrogen rather than gene
relative abundance. That is a FRESH-data / targeted-assay question, not a public-reanalysis one.

## Honest limits

- n≈34 healthy, WEIRD adults; capacity is **EC-level relative abundance** (compositional), summed over
  a small curated enzyme set — hypothesis-generating.
- Observational: diet is not randomized; diet and stool are paired by study-day, so `--lag` probes but
  cannot prove directionality.
- A within-person null does not refute between-person or interventional diet→function effects; it
  bounds what *daily, within-person, coarse-capacity* coupling looks like in this one cohort.
