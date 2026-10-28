# What we tried, what it bought us, and the tools that make it repeatable

A running record of the diet → gut → health investigation: every dataset and route attempted, the
honest outcome, and — the point of this note — the *value* each attempt produced even when the result
was null, plus the reusable tooling that lets the next question cost a fraction of the first.

The stance throughout is deliberate and worth stating, because it is what makes the nulls useful:

- **Public-data first.** Exhaust what can be reproduced from open data before asking anyone for
  controlled access or new sampling.
- **Reproducibility-first.** An effect only "counts" when it holds across independent cohorts (or,
  within a cohort, across people). A single-cohort signal is a hypothesis, not a result.
- **Honest nulls, bounded by positive controls.** A null is only informative if we can show the
  machinery *would* have detected a real effect. Every null below ships with that check.
- **Humble framing.** The contribution is rigor and reproducibility, not out-scaling a consortium.

---

## What we tried

| # | Dataset / route | What it tested | Outcome | What it bought us |
|---|---|---|---|---|
| 1 | **Italian StressDiet (n=113) + MELODY (n=47)**, reprocessed from raw 16S | diet → vaginal community state (CST-IV) | **Positive, replicated:** pooled animal-protein → CST-IV **OR 1.54 [1.08–2.20], p=0.017, I²=0%** | The one robust finding. MELODY alone was null (OR 1.30, p=0.44); pooling two cohorts recovered a consistent, significant signal — the method's first real payoff. |
| 2 | **cMD `westernized`** (Obregon-Tito 2015, n=58 stool shotgun) | does the only diet-proxy in curatedMetagenomicData track gut function? | glycogen & SCFA capacity *higher* in westernized (p<0.001), opposite the naive prior | Declared a **validated harness, not a finding** — and exposed the decisive structural fact: cMD's `westernized` never co-varies with a disease outcome in the same cohort, so cMD *alone* cannot do diet→function→health. Redirected the whole plan early. |
| 3 | **cMD function ↔ T2D** (Karlsson, Qin, Sankaranarayanan; 242 cases / 254 controls) | does gut SCFA / glycogen capacity associate with T2D, *consistently across cohorts*? | **Reproducible null:** SCFA pooled OR 1.19 [0.79–1.79], I²=61%; glycogen 1.05 — neither replicates | Showed the harness working as a guardrail: single cohorts would have mis-read (Sankaranarayanan alone SCFA OR 2.53, Qin 0.95); cross-cohort consistency collapsed it. Caught a per-study snapshot-date gotcha that now auto-resolves. |
| 4 | **SCFA panel broadening** (2 → 11 fermentation pathways) | does a richer functional definition rescue a signal? | Null held; heterogeneity *rose* (I² 49% → 61%) | Ruled out "we just defined the axis too narrowly." A negative that pre-empts an obvious reviewer objection. |
| 5 | **MetaCardis Zenodo archive** (record 5463864) | is the best gut cardiometabolic cohort an open mediation path? | Metabolome/microbiome feature tables are public (CC BY), but **FFQ diet + full clinical are request-only** | Confirmed it's the gut analog of PREDICT: features open, the diet+outcome join controlled. Saved us from assuming an open path that isn't. |
| 6 | **American Gut mediation** (HEI → α-diversity → BMI, n=839) | does gut diversity mediate diet's effect on BMI, on the one fully-open diet+microbiome+outcome set? | diet → BMI is **real** (total −0.715 [−1.058, −0.379]); **diversity mediates ~none** (NIE −0.036, CI includes 0, ~5%) | First end-to-end public mediation run. Established that the *signal is real but the mediator is too coarse* — the question that drove the next attempt. |
| 7 | **Johnson 2019 "Daily Sampling"** (34 adults × 17 days, daily ASA24 + daily shotgun; within-person) | does measured gut *function* (KEGG EC capacity), not diversity, track daily diet within a person? | **Reproducible null:** within-person r ≈ 0 for SCFA / estrobolome / glycogen; **I² ≈ 0** (not even personalized at this resolution) | The strongest test available — measured function, each person their own control — and it still washed out. Positive controls passed (within-person fiber↔kcal r=+0.58), so the null is real. Uniquely surfaced the estrobolome gene-level read (EC 3.2.1.31) that cMD couldn't give cheaply. |

---

## The value points — why this process is worth running

**1. Negative results with teeth.** Each null is bounded two ways: by a *positive control* (the machinery
demonstrably detects a real effect in the same data) and by *cross-unit consistency* (cohorts in #3,
people in #7). That turns "we didn't find anything" into "a coarse, public-data functional mediator does
not carry this signal, and here is the evidence it isn't a pipeline failure." That statement is reusable;
a bare null is not.

**2. The null is cheap because the tool is reusable.** Every attempt left behind a tested, documented
harness (synthetic demo → locked by a test → run on real data). The marginal cost of the next dataset is
a fetch/ingest script, not a new analysis. Seven routes share four core engines.

**3. Triangulation converges.** Three independent public substrates — cMD (cross-cohort function↔T2D),
American Gut (diversity mediation), Johnson 2019 (within-person daily function) — reach the *same*
conclusion from different angles. Agreement across independent designs is what lets us state the ceiling
with confidence rather than as one study's shrug.

**4. The process finds the decisive facts early.** The highest-leverage outputs were often structural, not
statistical: cMD's `westernized` never co-varies with outcome (#2); MetaCardis diet is request-only (#5);
cMD `gene_families` are 6M ID-only UniRef90 rows (#3→gene-level deferral). Each redirected effort before it
was wasted. Scoping the real metadata *is* a result.

**5. A bounded ceiling is a decision, not a dead end.** The cumulative finding —

> diet↔health signals are real (CST-IV OR 1.54; diet→BMI −0.715), **but the gut-function mediators
> buildable from public data — α-diversity *and* summed EC/pathway capacity — do not reproducibly carry
> them.**

— is precisely the justification for the one expensive move left: measured-function data with diet and
outcome in the same people (a controlled-access cohort, or fresh targeted sampling). The public work makes
that case concrete instead of assumed.

---

## The tools built (and how they let the research evolve)

### vmbdiet — the analysis engines (reused across all seven routes)

- **`associate.logistic_assoc`** — adjusted per-SD odds ratios (the field-standard diet↔microbiome model).
- **`associate.meta_analyze`** — inverse-variance meta-analysis (fixed / DerSimonian–Laird) with I²
  heterogeneity. The backbone of reproducibility-first: it pooled #1, flagged #3/#4 as non-replicating, and
  (via Fisher-z) pooled the per-person correlations in #7.
- **`functional`** — a curated, auditable axis dictionary (estrobolome / glycogen-starch / SCFA) grounded in
  MetaCyc pathways **and** EC/KO numbers, with `build_axis_feature_table`. Because it matches on EC, it ran
  unchanged on HUMAnN (cMD) *and* KEGG (Johnson 2019) tables — the same code, two data shapes.
- **`mediation`** — counterfactual NDE / NIE / proportion-mediated with cluster bootstrap (#6).
- **`maaslin`** — repeated-measures diet↔feature models for the clustered FRESH design (laptop engine +
  emitted R script).

### Reproducible-ingest pattern (one per dataset, data never redistributed)

`write_cmd_export_R` / `export_cmd_cohorts.R` (cMD), `fetch_americangut.py` (AGP),
`fetch_dietstudy.py` (Johnson 2019). Each turns a public source into a tidy, analysis-ready table; the raw
pulls are gitignored and reproduced on demand, so the repo stays light and redistribution-safe while every
result stays re-runnable.

### MiDAS — the layer that decides *whether* to pool before pooling

Cross-repo metadata indexing (MGnify, NMDC, ENA, dbGaP, cMD) plus a trust-aware reuse layer
(triangulation / provenance / comparability / self-validation). This is what scopes a question against the
*real* available metadata (as in #2 and #5) and gates which cohorts are comparable enough to meta-analyze —
turning "is there data for this?" into a checkable answer rather than a literature hunt.

### The loop that makes it an *evolving* process

```
question  →  scope the real metadata (MiDAS)  →  build a harness, lock it on synthetic data (a test)
          →  run on public data via a reproducible fetch  →  honest result (+ positive control)
          →  the harness persists  →  next question reuses it
```

The seventh route cost far less than the first because steps 3–5 were already built. That is the practical
meaning of "evolving research process": the questions change, the rigor is fixed, and the tooling
compounds.

---

## Honest limits of the exercise

- Public cohorts are observational, often self-reported, cross-sectional, and WEIRD-skewed; every result
  above is hypothesis-generating, and the one positive finding (#1) is an association, not causation.
- The functional axes are relative-abundance proxies summed over curated feature sets — deliberately coarse.
  A null at this resolution bounds the coarse story; it does not refute finer-resolution biology (specific
  butyrate producers, GUS structural subclasses, measured metabolites).
- "Public-data ceiling reached" is a claim about *these* open substrates for *this* mediation question, not
  a claim that no gut mechanism exists.
