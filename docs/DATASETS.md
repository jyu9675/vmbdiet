# Public datasets for diet ↔ vaginal microbiome

The hard requirement is **paired** data: per subject, both a vaginal microbiome profile *and* dietary
intake. Most vaginal-microbiome cohorts have no diet layer. These are the verified public exceptions
(accession status checked 2026-10; "verified live" = the archive record was opened directly).

## Primary targets (fully public, no controlled-access application)

| # | Cohort | Microbiome | Diet | N | Accession | Status |
|---|--------|-----------|------|---|-----------|--------|
| 1 | **Italian / Bologna** (Djusse/Zhu/Foschi 2025, *Front Cell Infect Microbiol*) | 16S amplicon | nutrient FFQ (supplement XLS) | 113 | **PRJNA1188525** (SRA) | ✅ verified live |
| 2 | **MELODY / IBD** (PLOS ONE 2024) | 16S V3–V4 | **HEI-2015** (3× 24-h recall, NDSR) | 48 | **PRJNA915128** (BioProject) | public per data-availability stmt |
| 3 | **MicrobeMom** (Microbiol Spectr 2024) | shotgun metagenomic | FFQ nutrient + amino acids | 118 | **PRJEB48251** (ENA) | ✅ verified live |

**Recommended first target: #1 (PRJNA1188525).** Preferred data type (16S), largest fully-open paired N,
non-pregnant reproductive-age (cleanest CST/Lactobacillus-dominance baseline), diet in the open Frontiers
supplement (`DataSheet1.xls` / `DataSheet2.xls`). **One thing to confirm first:** that the supplement
carries a subject-ID ↔ SRA-sample key to join diet to microbiome; if absent, email the corresponding author.
Article + supplements: https://www.frontiersin.org/articles/10.3389/fcimb.2025.1582283/full

- #2 is the direct **HEI-2015 → CST** cross-check (same metric as the Hawaii study), smaller + IBD confounding.
- #3 is a larger **replication** but shotgun (needs read→taxa processing to get CSTs) and a pregnancy cohort.

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
