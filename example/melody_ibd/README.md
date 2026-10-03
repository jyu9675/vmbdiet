# Second cohort — MELODY vaginal / IBD pregnancy cohort (PRJNA915128)

A second from-raw-reads run of `vmbdiet`, on an independent cohort with a different diet metric
(**HEI-2015**, the Healthy Eating Index — the same measure used by the Hawaii multi-ethnic study).
This is the honest companion to the Italian result: a smaller, pregnant cohort where the diet signals
**do not reach significance** — consistent with the original study's own conclusion.

## Source

| | |
|---|---|
| Paper | PLOS ONE 2024 — "Association of vaginal cytokines and dietary intake with IBD status and vaginal microbiota in pregnant individuals" — DOI [10.1371/journal.pone.0335178](https://doi.org/10.1371/journal.pone.0335178) (PMC12803450) |
| Cohort | 48 pregnant individuals, 3rd trimester — 23 with IBD (18 Crohn's, 5 UC) + 25 controls |
| Microbiome | vaginal 16S rRNA **V3–V4**, Illumina MiSeq 2×300, BioProject **PRJNA915128** (vaginal AMPLICON subset; the project also holds unrelated stool WGS runs, excluded) |
| Diet | 3× 24-h recall → **HEI-2015** + full nutrient panel |
| Join key | `sample_alias = <subject>_BSL1`; **per-sample diet + covariates are embedded in the BioSample metadata** (no separate supplement needed) — parsed in `scripts/fetch_data.sh` |

The deposited reads are already primer-trimmed, so the pipeline skips cutadapt (otherwise identical to the
Italian run). 1,072 ASVs → 128 taxa; 48/48 samples pass depth, 47 have complete diet.

## Reproduce

```bash
bash scripts/fetch_data.sh   # vaginal FASTQs + diet.tsv (HEI-2015 + nutrients from BioSample metadata)
bash scripts/pipeline.sh     # no primer trim (reads pre-trimmed) -> work/otutab.txt, zotus.fa, sintax.txt
python scripts/analyze.py    # CST + HEI-2015 / nutrient associations -> results/
```

## What we found (N=47, 12 CST-IV)

- **Community structure is realistic:** CST I 19 · III 14 · **IV 13** · II 1 · V 1 (37/48 *Lactobacillus*-dominant).
- **Diet → CST is null here.** Nothing is significant and nothing survives FDR:
  - HEI-2015 → CST-IV: OR 1.40 per SD (0.71–2.78), p=0.33
  - animal-protein density → CST-IV: OR 1.30 (0.67–2.51), p=0.44
- **But the animal-protein direction is concordant with the Italian cohort** (OR > 1 in both).
- **The alcohol → *Gardnerella* axis cannot be tested here:** this is a pregnancy cohort — alcohol intake is
  ≈ 0 in both groups (p=0.11), so there is no exposure contrast. No *Gardnerella* signal, as expected.

### Why a null is the right answer

This cohort is where power is weakest and the original authors themselves reported **no** significant
diet- or IBD-related difference in vaginal microbiota composition. Three structural reasons:
1. **N=47 with 12 CST-IV events** — a continuous per-SD OR of ~1.3 is undetectable at this size.
2. **Pregnancy flattens the exposure** — near-zero alcohol, and pregnancy strongly shifts the vaginal
   community toward *Lactobacillus* regardless of diet.
3. **IBD as an extra axis of variation** adds noise to a small sample.

Reporting this null alongside the positive Italian result is the point: the toolkit gives the same honest
answer the source study did, and the one directional signal that *is* comparable (animal protein) points the
same way in both cohorts.

See [`../../docs/italy-stressdiet-findings.md`](../../docs/italy-stressdiet-findings.md) for the companion
positive cohort and the cross-cohort summary.
