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

Pick a study (or pool several) that spans Westernized and Non-westernized samples. Then:

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
