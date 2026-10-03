# Public datasets for diet ↔ vaginal microbiome

The hard requirement is **paired** data: per subject, both a vaginal microbiome profile *and* dietary
intake. Most vaginal-microbiome cohorts have no diet layer. These are the verified public exceptions
(accession status checked 2026-10; "verified live" = the archive record was opened directly).

## Primary targets (fully public, no controlled-access application)

| # | Cohort | Microbiome | Diet | N | Accession | Status |
|---|--------|-----------|------|---|-----------|--------|
| 1 | **Italian / Bologna** (Foschi et al. 2025, *Front Cell Infect Microbiol*) | 16S amplicon | nutrient FFQ (supplement XLS) | 113 | **PRJNA1188525** (SRA) | ✅ **DONE** — [example/italy_stressdiet](../example/italy_stressdiet) |
| 2 | **MELODY / IBD pregnancy** (PLOS ONE 2024) | vaginal 16S V3–V4 | **HEI-2015** (3× 24-h recall) | 48 | **PRJNA915128** (BioProject) | ✅ **DONE** — [example/melody_ibd](../example/melody_ibd) |
| 3 | ~~**MicrobeMom** (Microbiol Spectr 2024)~~ | ~~shotgun~~ | FFQ | 118 | **PRJEB48251** (ENA) | ❌ **not usable — GUT/stool cohort, no vaginal samples** |

**Both runnable cohorts have now been analyzed** (see [italy-stressdiet-findings.md](italy-stressdiet-findings.md)).
Key practical lessons:
- **#1 (Italian)** — diet is in the open Frontiers supplement; `sample_alias = StressDiet_N` joins 1:1 to the
  supplement `sample` column. **Both published signals replicated** (animal protein → CST-IV; alcohol → *Gardnerella*).
- **#2 (MELODY)** — the diet (HEI-2015 + full nutrient panel) and covariates turned out to be embedded in the
  **BioSample metadata** (no supplement needed); `sample_alias = <subject>_BSL1`. Pull with `ena/browser/api/xml/<SAMN>`.
  Result: a **concordant-direction null** (underpowered pregnancy cohort; alcohol ≈ 0 so that axis is untestable).
- **#3 (MicrobeMom, PRJEB48251)** — verified at the ENA record to be **`human gut metagenome`** (1,451 runs, ~1.5 TB),
  a maternal-gut probiotic RCT with **no vaginal samples**. The earlier "vaginal replication" label was wrong; dropped.

## CST classification reference (VALENCIA)

France et al. 2020, *Microbiome* — the field-standard nearest-centroid CST classifier.
- Centroids (verified live): `https://raw.githubusercontent.com/ravel-lab/VALENCIA/master/CST_centroids_012920.csv`
  (header `sub_CST,Lactobacillus_iners,Lactobacillus_crispatus,Gardnerella_vaginalis,…`, ~245 taxa, sub-CSTs I-A…V).
- Classifier: `https://raw.githubusercontent.com/ravel-lab/VALENCIA/master/Valencia.py` · repo `ravel-lab/VALENCIA` (MIT).
- `vmbdiet.cst.assign_cst_valencia(abund, centroids)` scores against these centroids directly (Yue-Clayton), matching VALENCIA.

## Not currently usable

- **Hawaii multi-ethnic** (PMC11479099, HEI-2015 + 16S, N=40): data "available upon publication" but **no accession deposited yet** — revisit.
- **PIN prenatal** (PMC8881389, N≈634): no public accession (likely controlled).
- **Birmingham BV/nutrients** (PMC2663425, N=1521): **no sequencing at all** (Nugent score only) — useful only as a reference for which nutrients associate with BV.
- **MOMS-PI** (dbGaP phs001523): vaginal 16S but **no diet layer**.

## Run plan (real data)

1. **Diet:** download the cohort's nutrient table (supplement XLS for #1; derive HEI-2015 for #2).
2. **Microbiome → abundance table:** for 16S (#1, #2) process reads → ASV → species (DADA2/QIIME2) or use a provided ASV table; for shotgun (#3) map/profile → species. All laptop-feasible (16S is light).
3. **CST:** `vmbdiet cst` (dominant-taxon) or `assign_cst_valencia` with the centroids above.
4. **Associate:** `vmbdiet associate --outcome cst_IV|lacto_depleted --nutrients … --covariates age,race,bmi,…`.
5. Compare the adjusted ORs to the literature (fat↑ / animal-protein↑ / alcohol↑ risk; carbohydrate↑ / folate↑ / dairy↑ / fiber↑ protective).
