# Aim-3 functional contrast on curatedMetagenomicData (westernized proxy)

cMD doesn't carry detailed dietary intake, but it *is* strong on gut shotgun **function** (HUMAnN 3) and on
`westernized` — its best diet/lifestyle axis (non-westernized ≈ high-fiber, plant-rich). This harness tests
the Aim-3-flavored question **does gut functional capacity differ by westernization** (glycogen/starch and
SCFA biosynthesis — the fiber-fermentation expectation), using the same functional dictionary the gut–vagina
pipeline uses.

## One R step, then no more R

cMD's function tables ship via R, so export once (any machine with R):

```r
# BiocManager::install("curatedMetagenomicData"); vmbdiet.functional.write_cmd_export_R() emits this
library(curatedMetagenomicData); library(SummarizedExperiment)
se <- curatedMetagenomicData("2021-10-14.<Study>.pathway_abundance", dryrun=FALSE, counts=FALSE, rownames="long")[[1]]
write.table(as.data.frame(assay(se)),   "cmd_<Study>_pathway_abundance.tsv", sep="\t", quote=FALSE, col.names=NA)
write.table(as.data.frame(colData(se)), "cmd_<Study>_metadata.tsv",          sep="\t", quote=FALSE, col.names=NA)
```

Pick a study that spans both. Verified single-study contrasts (stool shotgun) in cMD:
**`Obregon-TitoAJ_2015`** (recommended — 22 Westernized US vs 36 Non-westernized Peru, n=58) or
**`RampelliS_2015`** (11 Italian vs 27 Hadza/Tanzania, n=38). The `westernized` column reads `Yes`/`No`.
Then:

```bash
python example/cmd_westernized/analyze.py \
    --pathways cmd_<Study>_pathway_abundance.tsv --meta cmd_<Study>_metadata.tsv
```

## Validate the harness offline (no R, no download)

```bash
python example/cmd_westernized/analyze.py      # runs the synthetic demo in make_demo.py
```

On the synthetic demo (planted non-westernized → higher fiber-fermentation capacity) the harness recovers
it: **glycogen_starch and SCFA capacity are significantly higher in Non-westernized samples** (Mann-Whitney
p < 0.001). That confirms the pipeline — HUMAnN pathway table → `build_axis_feature_table` → group contrast —
is correct; swap in a real cMD export to get the real result.

## Scope

- Pathway export covers **glycogen/starch** and **SCFA**. The **estrobolome** feature (β-glucuronidase) is
  gene-level (EC 3.2.1.31 / KO), so it needs the **gene-family** export (`data_type="gene_families"`,
  regrouped to EC/KO), not the pathway table.
- Westernization is a **proxy**, not measured diet — a strong, well-documented lifestyle axis, but coarser
  than the nutrient-intake exposures in the dedicated diet cohorts
  ([cross_cohort_meta](../cross_cohort_meta)). This is the functional/mechanism layer, complementary to the
  diet↔community association.

## Real-data run (Obregon-TitoAJ_2015) — result and conclusion

Run live on **`Obregon-TitoAJ_2015`** (stool shotgun, n=58: 22 US / westernized, 36 Peru /
non-westernized), via the cMD R export. Pathways matched: glycogen_starch = 3, scfa = 2.

| axis | mean (westernized, US) | mean (non-westernized, Peru) | p (MWU) | direction |
|---|---|---|---|---|
| glycogen_starch | 7.8e-4 | 5.8e-4 | <0.001 | higher in westernized |
| scfa | 1.4e-4 | 4.3e-5 | <0.001 | higher in westernized |

**Interpretation — treat as a validated harness, not a finding.** The contrast is statistically
significant, but it is a weak, narrow readout: tiny *relative* abundances summed over only 3 glycogen + 2
SCFA MetaCyc pathways (relative ≠ absolute capacity; compositional), on a coarse westernization proxy. The
direction is the *opposite* of the naive "non-westernized high-fiber → more SCFA" expectation — the
glycogen/starch direction is plausibly a high-starch Western diet, but the SCFA direction should not be
read as a real biological result from two butyrate pathways. This is deliberately **not** promoted to a
finding (unlike the clean, replicated diet→CST-IV result in [cross_cohort_meta](../cross_cohort_meta)); if
anything it is a cautionary illustration that pathway-relative-abundance proxies can contradict the fiber
narrative, which is exactly why the proposal needs real functional profiling + measured diet (the FRESH
cohort), not proxies.

**What this run *does* establish:** the end-to-end pipeline works on real cMD data (R export →
`build_axis_feature_table` → group contrast). A more defensible functional read, if ever pursued, is the
`gene_families` export (EC/KO) — a fuller SCFA capacity plus the estrobolome (β-glucuronidase) — but that
is out of scope here; we conclude the public-data functional exploration at the harness-validation stage.
