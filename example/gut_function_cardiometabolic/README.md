# Cross-cohort gut-function ↔ cardiometabolic reproducibility (Route A)

A reproducibility-first harness: does a gut functional **capacity** (SCFA / fiber fermentation,
glycogen/starch) associate with a cardiometabolic outcome *consistently across independent public
cohorts*? It computes an adjusted per-SD logistic association within each cohort, then pools the
per-cohort log-odds with `vmbdiet.meta_analyze` and flags an axis as **replicating** only if the
direction is the same in every cohort **and** the pooled estimate is significant. The contribution is
rigor — an effect only "counts" when it reproduces — not a new biological claim.

This is **Route A** of [`docs/diet-gut-health-landscape.md`](../../docs/diet-gut-health-landscape.md).

## Why this is an association, not a mediation

We scoped curatedMetagenomicData against its real `_sample.tsv` columns. cMD carries **no detailed
diet**, its only diet/lifestyle proxy is `westernized`, and — decisively — **`westernized` varies only
in healthy traditional-vs-Western cohorts, while every cMD cohort that has a cardiometabolic outcome
is uniformly `westernized=Yes`.** The exposure and the outcome never co-occur in one cohort, so cMD
*alone* cannot support a diet → function → health mediation. What it *can* support is this: function
↔ outcome reproducibility across cohorts. (BMI can't stand in as the diet exposure — it is downstream
of diet, an outcome/confounder.) The full diet→function→health mediation needs a cohort carrying
diet + function + outcome together — MetaCardis / American Gut / PREDICT — see the landscape doc.

## One R step (all cohorts), then no more R

cMD's function tables ship via R. [`export_cmd_cohorts.R`](export_cmd_cohorts.R) pulls all three T2D
cohorts in one pass — **`KarlssonFH_2013`** (n=145), **`QinJ_2012`** (n=363),
**`SankaranarayananK_2015`** (n=37) — writes each cohort's pathway table + metadata, and adds a
**coding-independent `t2d_case` 0/1 column** (cMD's disease coding varies by version, so the script
normalizes `study_condition`/`disease` itself) plus a `sex` alias. Its header lists the install
gotchas (Rtools45 for `bluster`; `rbiom` 1.0.3 downgrade; the ExperimentHub cache prompt). Run it on
any machine with R + `curatedMetagenomicData`:

```bash
Rscript export_cmd_cohorts.R
```

Then, in Python (no more R):

```bash
python analyze.py \
  --cohort Karlsson cmd_KarlssonFH_2013_pathway_abundance.tsv cmd_KarlssonFH_2013_metadata.tsv \
  --cohort Qin      cmd_QinJ_2012_pathway_abundance.tsv        cmd_QinJ_2012_metadata.tsv \
  --cohort Sankar   cmd_SankaranarayananK_2015_pathway_abundance.tsv cmd_SankaranarayananK_2015_metadata.tsv \
  --outcome t2d_case --covars age,sex,BMI
```

For an atherosclerotic-CVD axis instead, `JieZ_2017` (n=385) and `MetaCardis_2020_a` (n=1831) carry
`Atherosclerosis` / `Coronary Artery Disease` — point the script at those studies and map the case
label accordingly. (`analyze.py` also accepts a string outcome directly, e.g.
`--outcome disease --case "Type 2 Diabetes Mellitus" --control Healthy`, if you'd rather not use the
normalized column.)

## Validate the harness offline (no R, no download)

```bash
python analyze.py          # generates + runs the synthetic demo (make_demo.py)
```

The demo plants a modest, reproducible **SCFA-lower-in-T2D** signal across three cohorts and a null
glycogen axis. The harness recovers exactly that — SCFA pooled OR ≈ 0.46 [0.28–0.78], p ≈ 0.004,
replicates = **YES** (all three cohorts OR < 1); glycogen OR ≈ 0.85, p ≈ 0.20, replicates = no. Locked
by `tests/test_gut_cardiometabolic.py`. (The numbers above are synthetic — a pipeline check, not a
finding.)

## Scope and honest limits

- **Association, not mediation** (see above). Replication across independent cohorts controls for
  cohort-specific artefacts, but not for residual confounding beyond the adjusted covariates.
- Capacity is **relative-abundance** summed over a small MetaCyc pathway set (compositional); the axis
  dictionary is a curated starting set, not exhaustive. Treat results as hypothesis-generating.
- **Estrobolome** (β-glucuronidase) is gene-level (EC 3.2.1.31 / KO), so it needs the `gene_families`
  export, not the pathway table — it will not appear from pathway inputs.
- Pooling uses DerSimonian–Laird random effects by default; `I²` is reported so heterogeneity is
  visible. The MiDAS comparability/provenance layer can gate which cohorts are poolable before this
  step.
