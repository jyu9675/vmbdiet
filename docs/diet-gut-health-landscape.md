# Diet → gut → health: landscape + public-data scoping (2026-10-26)

A honest record of what we concluded when weighing a pivot from the gut–vagina axis toward
**how diet shapes the gut microbiome and, through it, cardiometabolic health** — and what public
data can actually support today. Written as a project record, not a result; nothing here is a
finding.

## 1. The question is real and under-answered — but the gap is methodological

Large consortia already map diet ↔ gut-microbiome ↔ health *associations* at scale (ZOE/PREDICT
ranked diet–microbe–health relationships across ~34k people). What is **not** established is the
*causal chain tested rigorously across independent cohorts*:

> diet → gut microbial **function** (SCFA / fiber fermentation, etc.) → cardiometabolic outcome,
> as a formal **mediation**, reproduced in more than one cohort.

Even the large studies stop at association and say causal inference needs prospective/interventional
designs; the 2024 community-scale metabolic-modelling work links *predicted* SCFA production to
cardiometabolic markers but is association-level, not mediation. **So the open niche is
methodological rigor** (formal mediation + cross-cohort reproducibility + confounder/batch control)
on public data — something an individual can contribute, rather than trying to out-scale a
consortium's proprietary cohort. This is exactly what `vmbdiet.mediation` + `vmbdiet.meta_analyze`
+ the MiDAS comparability/provenance layer are built to do.

## 2. The hard constraint: the data arm that's scarce is *detailed diet + health together*

Microbiomes are public almost everywhere. **Deeply-phenotyped diet + health data is usually
controlled or supplement-locked** — the same wall we hit on the vaginal axis, and it holds in
gut-health too:

- **PREDICT 1** (`AsnicarF_2021` in cMD) — microbiome/function **public** (1,098 stool samples;
  raw reads in SRA). Its detailed diet + cardiometabolic markers are **controlled** (released only
  to approved researchers). cMD carries only `bmi`, `disease=Healthy`, age, sex, `westernized=Yes`
  for it. → an *application target*, not an immediate open analysis.
- **MetaCardis** (`MetaCardis_2020_a` in cMD, n=1,831) — the best gut analog: microbiome/function
  **public** via cMD, a validated **FFQ** diet instrument + a large metabolic-marker panel + a full
  cardiometabolic-disease spectrum (CAD, T2D, heart failure, metabolic syndrome). But the detailed
  diet + markers are **not** in cMD's harmonized curation (which exposes only `disease`, `bmi`,
  `westernized`); they live in the source study, partly in a Zenodo pipeline archive
  ([zenodo 5463864](https://zenodo.org/record/5463864), contacts Forslund/Birkner — *check whether
  the archive's metadata includes FFQ-derived diet + markers or only features*), and 2025 linkage
  of French participants to the national health system (SNDS) means clinical outcomes are governed
  by French health-data rules. → the most promising lead, **access to be verified**.
- Genuinely-open, no application: **American Gut / AGP** (public VioScreen **FFQ** + 16S, large n —
  but mostly 16S, so function is inferred not measured) and small shotgun+diet cohorts (Johnson
  daily-sampling, German FFQ, Rocafort 24h-recall — already in MiDAS).

## 3. What curatedMetagenomicData can and cannot do (scoped against the real columns)

We read the actual `_sample.tsv` columns for the candidate cohorts. The curation schema's **only
diet/lifestyle field is `westernized`**, there are **no cardiometabolic-marker columns** anywhere
(cholesterol/glucose/HbA1c/BP/insulin…), and the only continuous phenotype co-present with function
across cohorts is **BMI**.

| cMD study | n | outcome(s) present | `westernized` | BMI | diet detail |
|---|---|---|---|---|---|
| KarlssonFH_2013 | 145 | T2D, glucose intolerance | Yes only | ✓ | none |
| QinJ_2012 | 363 | T2D | Yes only | ✓ | none |
| SankaranarayananK_2015 | 37 | T2D | Yes only | ✓ | none |
| JieZ_2017 | 385 | atherosclerosis (ACVD) | Yes only | ✓ | none |
| MetaCardis_2020_a | 1831 | CAD, T2D, HF, metabolic syndrome | Yes only | ✓ | none (in curation) |
| LeChatelierE_2013 | 292 | (obesity/BMI; all "Healthy") | Yes only | ✓ | none |
| YachidaS_2019 | 616 | CRC, adenoma | Yes only | ✓ | none |
| HMP_2019_ibdmdb / NielsenHB_2014 | 1627 / 396 | IBD | Yes only | ✓ | none |
| AsnicarF_2021 (PREDICT 1) | 1098 | all "Healthy" | Yes only | ✓ | none (in curation) |
| **RampelliS_2015** | 38 | all "Healthy" | **No/Yes** | — | none |
| **Obregon-TitoAJ_2015** | 58 | all "Healthy" | **No/Yes** | ✓ | none |

**The decisive finding:** `westernized` (cMD's only diet/lifestyle proxy) **varies only in healthy
traditional-vs-Western cohorts with no disease/cardiometabolic outcome**, and every cohort that
*has* a cardiometabolic outcome is uniformly `westernized=Yes`. The exposure and the outcome never
co-occur in the same cohort. **Therefore cMD alone cannot do the diet → function → health
mediation** — the diet arm is simply absent where the outcome exists. (Using BMI as the "diet" node
would be wrong — BMI is downstream of diet, an outcome/confounder, not the exposure.)

## 4. What's actually feasible, ranked

- **A — fully public, now (association, not mediation): cross-cohort FUNCTION ↔ cardiometabolic
  reproducibility.** Pool the T2D cohorts (Karlsson + Qin + Sankaranarayanan) and/or ACVD (Jie) +
  MetaCardis + the BMI cohorts, export HUMAnN function, and test whether specific functional
  capacities (SCFA, etc.) reproducibly associate with the outcome across *independent* cohorts,
  with the MiDAS comparability/provenance layer controlling batch and `meta_analyze` pooling. The
  contribution is the **reproducibility + trust rigor**, not a new association. This is the honest
  doable thing in cMD.
- **B — the true diet→function→health mediation: needs a cohort with all three.** Best lead
  **MetaCardis** (verify the Zenodo/source access to FFQ + markers); fully-open fallback **American
  Gut** (FFQ + 16S, function inferred); **PREDICT 1 / HPP** as parallel data-access applications.
- **C — the westernization → function contrast** (`example/cmd_westernized/`) is the only diet-proxy
  analysis cMD supports, and only in healthy cohorts — already run and concluded as a *validated
  harness, not a finding*.

**Bottom line:** the pivot is viable, but the win is methodological, and the first genuinely-public
result lives in route A (function↔outcome reproducibility) while route B's mediation waits on a
cohort that carries diet + function + outcome together.

Sources: PREDICT/ZOE and the MetaCardis dietary FFQ + data notes —
[MetaCardis FFQ validation (INSERM/HAL)](https://www.hal.inserm.fr/inserm-03970414),
[MetaCardis MICA catalogue](https://mica.mdc-berlin.de/study/metacardis),
[MetaCardis pipeline archive (Zenodo)](https://zenodo.org/record/5463864).
