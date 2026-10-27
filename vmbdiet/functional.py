"""
vmbdiet.functional — Aim-3 functional layer: HUMAnN 3 / eggNOG ingest + a mechanism-focused feature
dictionary for the gut-vagina axis.

Species-level association (Aims 1-2) says *who* is there; the axis hypothesis is about *what the gut
community can do* — specifically three mechanisms that could link diet to the vaginal niche:

  * estrobolome      — bacterial β-glucuronidase (GUS) deconjugates hepatic estrogens, raising
                       circulating estrogen -> vaginal epithelial glycogen -> lower pH favouring
                       L. crispatus (Takada 2023).
  * glycogen_starch  — glycogen / starch degradation capacity (the substrate side of the same story).
  * scfa             — short-chain fatty-acid (butyrate / acetate / propionate) biosynthesis.

This module (1) parses HUMAnN 3 output (gene families or pathway abundance, stratified or not) into a
clean feature table, and (2) subsets it to a curated, auditable dictionary of the identifiers that
actually appear in HUMAnN output — MetaCyc pathway IDs, level-4 EC numbers, KEGG orthologs, and name
patterns for raw UniRef tables. The matched table feeds straight into `vmbdiet.maaslin` (the same
random-effects models used for species), so Aim 3 is Aim 1-2 with a functional feature table swapped in.

Identifier provenance: EC numbers and MetaCyc pathway IDs are HUMAnN-native (regroup gene families with
`humann_regroup_table --groups uniref90_level4ec` / `uniref90_ko`). KOs listed here were checked against
KEGG (K01195 β-glucuronidase, K00925 acetate kinase, K00929 butyrate kinase). The dictionary is a
*starting* set to review with the PI, not an exhaustive one — see each category's `note`.

Important estrobolome caveat (Ervin et al. 2019): EC 3.2.1.31 / K01195 capture *total* bacterial GUS, but
only certain GUS structural classes (Loop 1, mini-Loop 1, FMN-binding) reactivate estrogens. KO/EC
annotation alone cannot isolate the estrogen-reactivating subset; that needs GUS sequence/structure
classification layered on top. Treat the estrobolome feature as "total GUS capacity", stated honestly.
"""
from __future__ import annotations
import re
import numpy as np
import pandas as pd

# Rows HUMAnN emits that are not real features.
_SPECIAL = {"UNMAPPED", "UNINTEGRATED", "UNGROUPED", "UNCLASSIFIED"}

# --------------------------------------------------------------------------- the dictionary

FUNCTIONAL_FEATURES: dict[str, dict] = {
    "estrobolome": {
        "label": "β-glucuronidase (estrogen deconjugation)",
        "pathways": [],
        "ecs": ["3.2.1.31"],
        "kos": ["K01195"],
        "cazy": ["GH2"],                      # reference only; HUMAnN does not regroup to CAZy natively
        "name_patterns": ["glucuronidase"],
        "note": "Total bacterial GUS capacity. KO/EC cannot isolate the estrogen-reactivating GUS "
                "loop-classes (Ervin 2019); refine with GUS structural classification for the "
                "estrobolome-specific signal.",
    },
    "glycogen_starch": {
        "label": "glycogen & starch degradation",
        "pathways": ["GLYCOCAT-PWY", "PWY-6731", "PWY-6737"],  # glycogen deg I; starch deg III, V
        "ecs": ["3.2.1.1", "3.2.1.3", "3.2.1.41", "2.4.1.25"],  # amylase, glucan-1,4-a-glucosidase, pullulanase, 4-a-glucanotransferase
        "kos": [],                             # rely on EC + pathway + name (amylase KO unverified -> omitted)
        "cazy": ["GH13"],                      # reference only
        "name_patterns": ["amylase", "glucosidase", "pullulanase", "glucanotransferase",
                           "glycogen phosphorylase"],
        "note": "Substrate side of the glycogen mechanism. GH13 is the α-amylase CAZy family; add CAZy "
                "regrouping (biobakery utility) if you want family-level resolution.",
    },
    "scfa": {
        "label": "short-chain fatty-acid biosynthesis (butyrate / acetate / propionate)",
        "pathways": [
            # butyrate
            "PWY-5676", "CENTFERM-PWY",           # acetyl-CoA->butanoate II; pyruvate->butanoate
            "PWY-5677", "P163-PWY", "PWY-6590",   # succinate->butanoate; Lys ferment->acetate+butanoate;
                                                  #   C. acetobutylicum acidogenic (butyrate/acetate)
            # propionate
            "P108-PWY", "PROPFERM-PWY", "PWY-5088",  # pyruvate->propanoate I; Ala->propanoate+acetate;
                                                     #   L-glutamate deg VIII -> propanoate
            # acetate / mixed-acid
            "PWY-5100", "PWY-7254", "FERMENTATION-PWY",  # pyruvate->acetate+lactate II; TCA VII
                                                         #   acetate-producers; mixed-acid fermentation
        ],
        "ecs": ["2.8.3.8", "2.7.2.7", "2.7.2.1", "2.3.1.8", "5.4.99.2", "2.8.3.1"],
        # butyryl-CoA:acetate CoA-transferase, butyrate kinase, acetate kinase, phosphate
        # acetyltransferase, methylmalonyl-CoA mutase (propionate), propionate CoA-transferase
        "kos": ["K00925", "K00929"],           # acetate kinase (ackA), butyrate kinase (buk)
        "cazy": [],
        "name_patterns": ["butyrate kinase", "butyryl-coa", "acetate kinase",
                          "phosphate acetyltransferase", "methylmalonyl-coa mutase", "propionate"],
        "note": "Terminal SCFA enzymes + a broadened bacterial-fermentation pathway panel: butyrate "
                "(acetyl-CoA / pyruvate / succinate / lysine routes), propionate (pyruvate / alanine / "
                "glutamate routes) and acetate / mixed-acid. Excludes solventogenic (butanol/acetone/"
                "isopropanol), lactate-only, methanogenic and aromatic-degradation pathways, which are not "
                "SCFA biosynthesis. EC 5.4.99.2 / 2.8.3.1 remain the robust propionate enzyme markers.",
    },
}

CATEGORIES = tuple(FUNCTIONAL_FEATURES.keys())


# --------------------------------------------------------------------------- HUMAnN parsing

def feature_id(raw: str) -> str:
    """Bare feature id from a HUMAnN row label: strips the `|g__...s__...` stratification and the
    `: description` suffix. 'PWY-5676: acetyl-CoA...|g__X.s__Y' -> 'PWY-5676'."""
    base = str(raw).split("|", 1)[0]
    base = base.split(":", 1)[0]
    return base.strip()


def read_humann_table(path_or_df, stratified: str = "total", drop_special: bool = True) -> pd.DataFrame:
    """Load a HUMAnN 3 gene-family or pathway-abundance table (single-sample or `humann_join_tables`
    merged). Rows are features, columns are samples.

    stratified: 'total' keeps community totals only (rows without `|`, the usual association unit),
                'species' keeps only per-species rows, 'keep' keeps everything.
    drop_special: remove UNMAPPED / UNINTEGRATED / UNGROUPED / UNCLASSIFIED.
    Adds a `feature_id` column (bare id) and keeps the original label as the index.
    """
    if isinstance(path_or_df, pd.DataFrame):
        df = path_or_df.copy()
    else:
        df = pd.read_csv(path_or_df, sep="\t", header=0, index_col=0)
    df.index = df.index.map(str)
    strat = df.index.str.contains(r"\|", regex=True)
    if stratified == "total":
        df = df[~strat]
    elif stratified == "species":
        df = df[strat]
    ids = pd.Index([feature_id(x) for x in df.index])
    if drop_special:
        df = df[~ids.isin(_SPECIAL)]
    df.insert(0, "feature_id", [feature_id(x) for x in df.index])
    return df


def detect_level(ids) -> str:
    """Guess the identifier level of a set of feature ids: 'ko', 'ec', 'pathway', or 'uniref'."""
    ids = [str(i) for i in ids]
    n = max(len(ids), 1)
    ko = sum(bool(re.fullmatch(r"K\d{5}", i)) for i in ids)
    ec = sum(bool(re.fullmatch(r"\d+\.\d+\.\d+\.\d+", i)) for i in ids)
    uni = sum(i.startswith("UniRef") for i in ids)
    pwy = sum(("PWY" in i) or i.endswith("-PWY") for i in ids)
    best = max([("ko", ko), ("ec", ec), ("uniref", uni), ("pathway", pwy)], key=lambda t: t[1])
    return best[0] if best[1] > 0.1 * n else "uniref"


# --------------------------------------------------------------------------- matching

def _dict_ids(level: str, categories):
    key = {"ko": "kos", "ec": "ecs", "pathway": "pathways"}.get(level)
    out = {}
    for cat in categories:
        spec = FUNCTIONAL_FEATURES[cat]
        out[cat] = set(spec.get(key, []) or []) if key else set()
    return out


def match_features(table: pd.DataFrame, categories=None, level: str = "auto") -> pd.DataFrame:
    """Label the rows of a parsed HUMAnN table that belong to each axis category.

    Matches by exact feature-id for pathway/EC/KO tables; for raw UniRef tables (no EC/KO), matches the
    row's name against each category's `name_patterns`. Returns a long frame: feature, feature_id,
    category, matched_by, + the original abundance columns.
    """
    categories = list(categories or CATEGORIES)
    if "feature_id" not in table.columns:
        table = read_humann_table(table)
    if level == "auto":
        level = detect_level(table["feature_id"])
    id_sets = _dict_ids(level, categories)
    abund_cols = [c for c in table.columns if c != "feature_id"]
    rows = []
    for label, r in table.iterrows():
        fid = r["feature_id"]
        name = str(label).lower()
        for cat in categories:
            matched_by = None
            if fid in id_sets[cat]:
                matched_by = level
            elif any(pat.lower() in name for pat in FUNCTIONAL_FEATURES[cat]["name_patterns"]):
                matched_by = "name"
            if matched_by:
                rows.append({"feature": label, "feature_id": fid, "category": cat,
                             "matched_by": matched_by, **{c: r[c] for c in abund_cols}})
    return pd.DataFrame(rows)


def build_axis_feature_table(table, categories=None, level: str = "auto",
                             aggregate: str = "sum") -> pd.DataFrame:
    """Turn a HUMAnN table into a **samples x features** table ready for `vmbdiet.maaslin`.

    aggregate='sum'  -> one column per category (summed capacity, the usual axis score);
    aggregate='none' -> one column per matched feature, named '<category>|<feature_id>'.
    HUMAnN tables are features x samples, so the result is transposed to samples x features.
    """
    if not isinstance(table, pd.DataFrame) or "feature_id" not in getattr(table, "columns", []):
        table = read_humann_table(table)
    m = match_features(table, categories=categories, level=level)
    abund_cols = [c for c in table.columns if c != "feature_id"]
    if m.empty:
        return pd.DataFrame(index=abund_cols)
    if aggregate == "sum":
        wide = m.groupby("category")[abund_cols].sum().T           # samples x category
    else:
        m = m.copy()
        m["col"] = m["category"] + "|" + m["feature_id"].astype(str)
        wide = m.groupby("col")[abund_cols].sum().T               # samples x feature
    wide.index.name = "sample"
    return wide


# --------------------------------------------------------------------------- curatedMetagenomicData export (R)

def write_cmd_export_R(out_path: str, study: str, data_type: str = "pathway_abundance",
                       date: str = "2021-10-14") -> str:
    """Emit an R snippet that pulls a curatedMetagenomicData study's functional table + sample metadata
    to TSVs this module can then read — the public-data prototype path for Aim 3 (cMD ships pre-computed
    HUMAnN 3 pathway / gene-family tables, so no read-level compute is needed). Run on a machine with R +
    Bioconductor `curatedMetagenomicData`."""
    script = f'''#!/usr/bin/env Rscript
# Auto-generated by vmbdiet.functional — export a curatedMetagenomicData study's HUMAnN 3 output.
# install once: BiocManager::install("curatedMetagenomicData")
suppressMessages(library(curatedMetagenomicData))
suppressMessages(library(SummarizedExperiment))

se <- curatedMetagenomicData("{date}.{study}.{data_type}", dryrun = FALSE, counts = FALSE,
                             rownames = "long")[[1]]
feat <- as.data.frame(assay(se))                     # features x samples (HUMAnN layout)
meta <- as.data.frame(colData(se))                   # samples x metadata (diet cols, covariates)

write.table(feat, "cmd_{study}_{data_type}.tsv", sep = "\\t", quote = FALSE, col.names = NA)
write.table(meta, "cmd_{study}_metadata.tsv",     sep = "\\t", quote = FALSE, col.names = NA)
cat("wrote cmd_{study}_{data_type}.tsv and cmd_{study}_metadata.tsv\\n")
# then: vmbdiet functional --humann cmd_{study}_{data_type}.tsv ...  ->  vmbdiet maaslin ...
'''
    with open(out_path, "w", newline="\n", encoding="utf-8") as fh:
        fh.write(script)
    return out_path
