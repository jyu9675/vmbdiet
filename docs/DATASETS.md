# Public datasets for diet ↔ vaginal microbiome

The hard requirement is **paired** data: per subject, both a vaginal microbiome profile *and* dietary
intake. Most vaginal-microbiome cohorts have no diet layer. These are the verified public exceptions
(accession status checked 2026-10; "verified live" = the archive record was opened directly).

## Primary targets (fully public, no controlled-access application)

| # | Cohort | Microbiome | Diet | N | Accession | Status |
|---|--------|-----------|------|---|-----------|--------|
| 1 | **Italian / Bologna** (Foschi et al. 2025, *Front Cell Infect Microbiol*) | 16S amplicon | nutrient FFQ (supplement XLS) | 113 | **PRJNA1188525** (SRA) | ✅ **DONE** — [example/italy_stressdiet](../example/italy_stressdiet) |
| 2 | **MELODY / IBD pregnancy** (PLOS ONE 2024) | vaginal 16S V3–V4 | **HEI-2015** (BioSample metadata) | 48 | **PRJNA915128** (BioProject) | ✅ **DONE** — [example/melody_ibd](../example/melody_ibd) |
| 3 | **IATA artificial-insemination** (Food & Function 2026, DOI 10.1039/d5fo04208a) | vaginal 16S amplicon | **MEDI-LITE** (Mediterranean) | 106 | **PRJNA1234600** (SRA) | 🟡 reads **staged**; diet table behind RSC paywall |
| 4 | **MicrobeMom** (ASM Spectrum 2024, PMC11537119) | **283 vaginal** 16S (+ stool/oral/milk) | 3-day diary → macro/micro + amino acids | ~118 | **PRJEB48251** (ENA) | 🟡 vaginal reads public; diet **on request** |

**Cohorts #1 and #2 analyzed** (see [italy-stressdiet-findings.md](italy-stressdiet-findings.md)); #3 is the top open lead; #4 corrects an earlier error.
Key practical lessons:
- **#1 (Italian)** — diet in the open Frontiers supplement; `sample_alias = StressDiet_N` joins 1:1. **Both published signals replicated** (animal protein → CST-IV; alcohol → *Gardnerella*).
- **#2 (MELODY)** — diet (HEI-2015 + full nutrient panel) + covariates embedded in the **BioSample metadata** (`ena/browser/api/xml/<SAMN>`); `sample_alias = <subject>_BSL1`. Result: a **concordant-direction null**.
- **#3 (IATA, PRJNA1234600)** — 106 vaginal 16S samples public & verified; study is explicitly Mediterranean-diet (MEDI-LITE) × vaginal microbiota × AI pregnancy success. Reads downloaded/staged in this repo's workflow. **Blocker:** the per-swab MEDI-LITE table lives only in the RSC Electronic Supplementary Information (`d5fo04208a`); the article is not open-access and RSC returns 403 to automated fetches. **To finish:** obtain the ESI (institutional access or author request), map it to `sample_alias`/`source_material_id` ("swab N"), then run `vmbdiet associate`.
- **#4 (MicrobeMom, PRJEB48251)** — **correction of an earlier note.** PRJEB48251 is a mixed-body-site cohort that **does contain 283 `Vaginal_microbiome` samples** (aliases `PV###_{E/L/1M}`; the earlier "no vaginal samples / gut-only" note was wrong — it mis-read `scientific_name` and missed `environment_material`). The real blocker is the diet join: the 3-day-diary nutrients are **not** in a per-subject supplement table and the data-availability statement routes them to the corresponding author. Runnable if the author supplies a PV-ID-keyed nutrient table.

## CST classification reference (VALENCIA)

France et al. 2020, *Microbiome* — the field-standard nearest-centroid CST classifier.
- Centroids (verified live): `https://raw.githubusercontent.com/ravel-lab/VALENCIA/master/CST_centroids_012920.csv`
  (header `sub_CST,Lactobacillus_iners,Lactobacillus_crispatus,Gardnerella_vaginalis,…`, ~245 taxa, sub-CSTs I-A…V).
- Classifier: `https://raw.githubusercontent.com/ravel-lab/VALENCIA/master/Valencia.py` · repo `ravel-lab/VALENCIA` (MIT).
- `vmbdiet.cst.assign_cst_valencia(abund, centroids)` scores against these centroids directly (Yue-Clayton), matching VALENCIA.

## Controlled-access / request-only (strong vaginal+diet cohorts behind a gate)

These have real vaginal microbiome **and** diet data but require a data application or author request; not runnable now, listed so the access path is on record (verified 2026-10).

- **Human Phenotype Project / Weizmann 10K** — the richest diet×vaginal pairing that exists: gut + **vaginal** + oral microbiome with app-based diet logs, FFQ, and CGM (~10k with diet + shotgun; Nat Med 2026). ENA projects **PRJEB85771 / PRJEB85945** verify but return **zero public runs** (managed access). Both reads and diet are gated. Apply: humanphenotypeproject.org/data-access (info@pheno.ai), bona-fide-researcher DAA.
- **Emory African-American Microbiome in Pregnancy** (Dunlop/Corwin/Brennan). Vaginal 16S reads **public**: **PRJNA725416** (436 `human vaginal metagenome`); **PRJNA553594** is the *oral* arm, not vaginal. Phenotypes (whether diet is included is **unverified** — protocol emphasizes psychosocial measures) are in a **dbGaP** controlled study (exact `phs` not resolved; it is **not** phs001523, which is MOMS-PI/VCU). Apply via dbGaP.
- **PIN vaginal / "Race & the Vaginal Microbiome & sPTB"** (mSystems 2022, PMC9238383). **PRJNA694098 — 824 vaginal 16S, public reads.** Nested in the Pregnancy-Infection-and-Nutrition (PIN) cohort, which holds FFQ/nutrition, but diet is **not** in the public BioSample — request from PIN investigators (UNC CPC).
- **Isala** (Belgium citizen-science; Nat Microbiol 2023). Reads public: ENA **PRJEB50407** (3,345 vaginal 16S). Per-subject metadata **including the diet FFQ is EGA-controlled: EGAD00001009890** (DAA + DAC, ~2–3 mo).
- **HCL / Tuddenham–Brotman** (PMC7387193 2020; PMC6806504 2019) — vaginal 16S V3–V4, Block Brief 2000 FFQ (macronutrient-fiber, betaine), n≈104–121, but data-availability says reads "**will be** released" (no accession found) — email authors / re-check SRA.

## Not usable (fail a hard requirement)

- **Hawaii multi-ethnic** (PMC11479099, HEI-2015 + full-length V2–V9 16S, N=40): **no accession deposited yet** ("following acceptance") — revisit for a PRJNA/GEO ID; pipeline needs full-length (no pair-merge) handling.
- **Birmingham BV/nutrients** (PMC2663425, N=1521): **no sequencing** (Nugent score only) — reference for which nutrients associate with BV, not a reanalysis target.
- **MOMS-PI** (dbGaP phs001523, VCU): large vaginal 16S but **no diet layer**.
- **mSphere daily-fluctuations** (PRJNA637322, 1,101 vaginal 16S, public): diet is only vegetarian/non-veg binary and "on request" — too coarse.
- **IMPACT BCN MedDiet RCT** (AJCN 2025, PMC12674037): excellent diet (preg-MEDAS + 151-item FFQ) + vaginal 16S N≈351, but **no accession issued yet** (placeholder in data-availability) + diet behind a DAA — re-check for the BioProject.
- Chinese cohort (PMC11836416, 6,755, 16S) — no dietary intake. Iranian BV-diet papers / NHANES BV-diet — FFQ+BV but **no sequencing**.

## Run plan (real data)

1. **Diet:** download the cohort's nutrient table (supplement XLS for #1; derive HEI-2015 for #2).
2. **Microbiome → abundance table:** for 16S (#1, #2) process reads → ASV → species (DADA2/QIIME2) or use a provided ASV table; for shotgun (#3) map/profile → species. All laptop-feasible (16S is light).
3. **CST:** `vmbdiet cst` (dominant-taxon) or `assign_cst_valencia` with the centroids above.
4. **Associate:** `vmbdiet associate --outcome cst_IV|lacto_depleted --nutrients … --covariates age,race,bmi,…`.
5. Compare the adjusted ORs to the literature (fat↑ / animal-protein↑ / alcohol↑ risk; carbohydrate↑ / folate↑ / dairy↑ / fiber↑ protective).
