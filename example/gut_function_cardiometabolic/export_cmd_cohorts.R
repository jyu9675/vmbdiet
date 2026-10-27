#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------------------------------
# Export the three cMD T2D cohorts' HUMAnN 3 pathway tables + normalized metadata to TSVs that
# example/gut_function_cardiometabolic/analyze.py reads. Run on a machine with R + Bioconductor's
# curatedMetagenomicData. This is the ONLY R step; everything downstream is Python.
#
# INSTALL (R 4.5 + Bioconductor 3.22, Windows — hard-won, see vmbdiet docs):
#   install.packages("BiocManager")
#   # Rtools45 must be installed (compiler for `bluster`) — that's NOT RStudio; it's the toolchain:
#   #   https://cran.r-project.org/bin/windows/Rtools/rtools45/rtools.html
#   BiocManager::install("curatedMetagenomicData", update = FALSE, ask = FALSE)
#   # rbiom 2.x dropped `unifrac` which mia/cMD imports -> downgrade:
#   install.packages("remotes"); remotes::install_version("rbiom", "1.0.3")
#   # The first curatedMetagenomicData() call prompts to create the ExperimentHub cache dir -> say YES
#   # (don't paste other commands into that prompt). Re-run this script after the cache exists.
# ---------------------------------------------------------------------------------------------------
suppressMessages(library(curatedMetagenomicData))
suppressMessages(library(SummarizedExperiment))

studies <- c("KarlssonFH_2013", "QinJ_2012", "SankaranarayananK_2015")  # stool shotgun, T2D vs healthy

for (study in studies) {
  cat("\n==== ", study, " ====\n", sep = "")
  # Resolve the study's pathway_abundance resource across cMD snapshots -- studies differ by date
  # (Karlsson/Qin are 2021-10-14 but Sankaranarayanan is only 2021-03-31), so take the LATEST match
  # rather than hardcode a date.
  avail <- curatedMetagenomicData(paste0("\\.", study, "\\.pathway_abundance$"), dryrun = TRUE)
  if (length(avail) < 1) { cat("  SKIP: no pathway_abundance resource for ", study, "\n", sep = ""); next }
  res <- tail(sort(avail), 1)
  cat("  resource: ", res, "\n", sep = "")
  se  <- curatedMetagenomicData(res, dryrun = FALSE, counts = FALSE, rownames = "long")[[1]]

  feat <- as.data.frame(assay(se))        # features x samples (HUMAnN layout) — relative abundance
  md   <- as.data.frame(colData(se))      # samples x metadata

  # --- normalize the outcome to a coding-independent 0/1 column (cMD coding varies by version) ---
  cond <- if ("study_condition" %in% colnames(md)) as.character(md$study_condition) else rep(NA, nrow(md))
  dis  <- if ("disease"         %in% colnames(md)) as.character(md$disease)         else rep(NA, nrow(md))
  tag  <- paste(cond, dis)
  is_case <- grepl("T2D|type 2 diab", tag, ignore.case = TRUE)
  is_ctrl <- grepl("control|healthy", tag, ignore.case = TRUE) & !is_case
  md$t2d_case <- ifelse(is_case, 1L, ifelse(is_ctrl, 0L, NA_integer_))
  if ("gender" %in% colnames(md)) md$sex <- md$gender        # analyze.py's default covars use `sex`

  write.table(feat, paste0("cmd_", study, "_pathway_abundance.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)
  write.table(md,   paste0("cmd_", study, "_metadata.tsv"),
              sep = "\t", quote = FALSE, col.names = NA)

  cat("  wrote cmd_", study, "_pathway_abundance.tsv (", nrow(feat), " features x ", ncol(feat),
      " samples)\n", sep = "")
  cat("  T2D cases =", sum(md$t2d_case == 1, na.rm = TRUE),
      " controls =", sum(md$t2d_case == 0, na.rm = TRUE),
      " unclassified =", sum(is.na(md$t2d_case)), "\n")
  cat("  raw study_condition / disease coding seen:\n")
  print(sort(unique(tag)))
}

cat("\nDone. Next (Python, no more R):\n",
    "  python analyze.py \\\n",
    "    --cohort Karlsson cmd_KarlssonFH_2013_pathway_abundance.tsv cmd_KarlssonFH_2013_metadata.tsv \\\n",
    "    --cohort Qin      cmd_QinJ_2012_pathway_abundance.tsv        cmd_QinJ_2012_metadata.tsv \\\n",
    "    --cohort Sankar   cmd_SankaranarayananK_2015_pathway_abundance.tsv cmd_SankaranarayananK_2015_metadata.tsv \\\n",
    "    --outcome t2d_case --covars age,sex,BMI\n", sep = "")
