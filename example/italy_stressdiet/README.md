# Real-data validation — Italian "StressDiet" cohort (PRJNA1188525)

An end-to-end, from-raw-reads reproduction of a published diet ↔ vaginal-microbiome study,
run entirely with **vmbdiet** plus a lightweight open-source 16S pipeline. Nothing here reuses the
authors' processed tables: the ASV table is rebuilt from the raw SRA FASTQs, CSTs are assigned by
vmbdiet, and the diet associations are computed by `vmbdiet.associate`.

## Source

| | |
|---|---|
| Paper | Laghi, Foschi, et al. *Front. Cell. Infect. Microbiol.* 2025 — DOI [10.3389/fcimb.2025.1582283](https://doi.org/10.3389/fcimb.2025.1582283) (PMC12213695) |
| Cohort | 113 reproductive-age women, Bologna, Italy; vaginal swabs + diet + ¹H-NMR metabolome |
| Microbiome | 16S rRNA **V3–V4** (341F/805R), Illumina MiSeq 2×300, BioProject **PRJNA1188525** |
| Diet | semi-quantitative FFQ (EPIC), 17 nutrient variables + MEDI-LITE Mediterranean score |
| Join key | SRA `sample_alias = StressDiet_N` ↔ supplement `sample = N` (perfect 1:1, N=113) |
| License | article + supplement CC-BY 4.0 (diet/metabolite tables redistributed here with attribution) |

## Pipeline (all open-source, laptop-feasible)

```
raw FASTQ ─cutadapt─▶ primer-trimmed ─vsearch─▶ merge(84–89%) ─▶ EE<1.0 filter
         ─▶ UNOISE3 ASVs (+UCHIME3) ─▶ OTU table ─▶ SINTAX/RDP genus
         ─▶ Lactobacillus species by usearch_global vs curated RefSeq 16S (Tier-2)
         ─▶ vmbdiet: CST (dominant-taxon) ─▶ logistic / Spearman / PERMANOVA
```

- primer trim: **cutadapt 5.2**; merge/denoise/taxonomy: **vsearch 2.32.0**
- genus reference: RDP SINTAX `rdp_16s_v16`; species reference: `data/vaginal_species_db.fa`
  (NCBI RefSeq `NR_` 16S for *L. crispatus/iners/gasseri/jensenii*, *G. vaginalis*, *F. vaginae*, + others).
  The 4 CST-defining *Lactobacillus* species are 94–97% identical across V3–V4 — enough margin for
  confident per-ASV assignment.
- 5.1M merged reads → **1,312 ASVs → 163 taxa**; every sample kept (depth ≥ 1000).

## Reproduce

```bash
# deps (once): cutadapt, and the vsearch 2.32 Windows/Linux binary on PATH or under tools/;
# RDP SINTAX reference at ref/rdp_16s_v16.fa.gz
bash scripts/fetch_data.sh     # supplement (nutrients/metabolites) + ENA manifest + ~3 GB FASTQs
bash scripts/pipeline.sh       # FASTQ -> work/otutab.txt, work/zotus.fa, work/sintax.txt
python scripts/analyze.py      # -> results/*.tsv, results/summary.json  (CST + associations)
```

## What we found (N=113)

**Community structure (sanity check — matches a reproductive-age cohort):**
CST **I** (*L. crispatus*) 49 · **III** (*L. iners*) 22 · **IV** (diverse/BV) 35 · **V** 4 · **II** 3.
Mean abundances: *L. crispatus* 41%, *L. iners* 22%, *G. vaginalis* 6.3% (81% prevalence) — textbook.

**Both of the paper's headline diet signals independently replicate:**

| Finding | Our result (independent reprocessing) | Direction vs paper |
|---|---|---|
| **Animal protein → CST-IV/BV** | energy-adjusted density **OR 1.66 per SD (1.09–2.53), p=0.019** (Mann-Whitney on density p=0.026) | ✅ replicates |
| **Alcohol → *Gardnerella*** | Spearman **rho +0.30, p=0.001** (q=0.068); ***Leptotrichia* rho +0.36, q=0.012** survives FDR | ✅ replicates + extends |

**Methodological nuance (and why it matters):** on *raw* intake, vegetable protein and fiber look
protective (lower in CST-IV, p=0.005 / 0.021) — but total energy is also lower in CST-IV (p=0.016),
so that signal is a quantity artifact. Switching to **energy density (per 1000 kcal)** removes it and
surfaces the real diet-*quality* effect: a higher **animal-protein fraction** tracks with CST-IV.
This is exactly the confounder energy adjustment exists to catch.

**Honest limits:** N=113 / 35 CST-IV events (modest power); the animal-protein association is
nominally significant but q=0.14 after BH-FDR across 7 nutrient densities; CSTs are cross-sectional.
PERMANOVA on community ~ CST-IV is strong (pseudo-F≈22, p=0.001) but partly definitional.

See [`../../docs/italy-stressdiet-findings.md`](../../docs/italy-stressdiet-findings.md) for the full write-up.
