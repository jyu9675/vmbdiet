# Functional layer (Aim 3)

`vmbdiet.functional` turns HUMAnN 3 output into a mechanism-focused feature table for the gut–vagina
axis, then hands it to [`vmbdiet.maaslin`](maaslin-bridge.md) — so Aim 3 (diet ↔ microbial *function*)
reuses the Aim 1–2 association engine with a functional feature table swapped in for the species table.

## Why these three mechanisms

The axis hypothesis (Takada 2023) is that the gut community shapes the vaginal niche through specific
functions, not just taxonomy:

| category | mechanism | why it matters to the vagina |
|---|---|---|
| **estrobolome** | bacterial β-glucuronidase (GUS) deconjugates hepatic estrogens | ↑ circulating estrogen → ↑ vaginal epithelial glycogen → ↓ pH → favours *L. crispatus* |
| **glycogen_starch** | glycogen / starch degradation | the substrate side of the same glycogen story |
| **scfa** | butyrate / acetate / propionate biosynthesis | gut SCFAs modulate mucosal immunity and the gut environment |

## The dictionary

`FUNCTIONAL_FEATURES` maps each category to the identifiers that actually appear in HUMAnN 3 output, at
every level HUMAnN can emit:

- **MetaCyc pathway IDs** (pathway-abundance tables): `GLYCOCAT-PWY` (glycogen deg I), `PWY-6731`/`PWY-6737`
  (starch deg III/V), `PWY-5676`/`CENTFERM-PWY` (butyrate).
- **Level-4 EC numbers** (`humann_regroup_table --groups uniref90_level4ec`): β-glucuronidase `3.2.1.31`;
  amylase/glucosidase/pullulanase/glucanotransferase `3.2.1.1`/`3.2.1.3`/`3.2.1.41`/`2.4.1.25`; SCFA
  `2.8.3.8`/`2.7.2.7`/`2.7.2.1`/`2.3.1.8` and propionate `5.4.99.2`/`2.8.3.1`.
- **KEGG orthologs** (`--groups uniref90_ko`), verified against KEGG: `K01195` (GUS), `K00925` (acetate
  kinase), `K00929` (butyrate kinase).
- **name patterns** for raw UniRef tables that have not been regrouped (`glucuronidase`, `amylase`,
  `acetate kinase`, …).
- **CAZy families** (GH2 for GUS, GH13 for amylases) are listed for reference; HUMAnN does not regroup to
  CAZy natively.

`match_features()` auto-detects the table's level and matches accordingly; `build_axis_feature_table()`
sums each category to a capacity score (or keeps per-feature columns) as a **samples × features** table
for the association step.

**This dictionary is a starting, auditable set — review it with the PI, not a finished ontology.** Every
ID's provenance is in the module docstring and each category carries a `note`.

## The estrobolome caveat (important)

EC `3.2.1.31` / KO `K01195` capture **total** bacterial β-glucuronidase. But only certain GUS structural
classes (Loop 1, mini-Loop 1, FMN-binding) reactivate estrogens (Ervin et al. 2019, *J Biol Chem*);
KO/EC annotation alone cannot isolate that subset. So the estrobolome feature means "total GUS capacity".
Refining it to the estrogen-reactivating GUS requires sequence/structure classification layered on the
UniRef hits — a natural extension once the marginal signal is established.

## Normalization

A **full** functional profile (hundreds of ECs/pathways) is compositional → feed it to `maaslin` with the
default TSS, then subset results to the axis features. A **pre-aggregated axis-score** table (a few
capacity sums) is *not* compositional → use `--normalize none` (absolute, log-transformed), or TSS will
close the few columns and force spurious negative associations. The functional demo uses `none`.

## Public prototype path

`write_cmd_export_R()` emits an R snippet that pulls a **curatedMetagenomicData** study's pre-computed
HUMAnN 3 pathway/gene table + sample metadata to TSV — no read-level compute, nothing heavier than a
download. `vmbdiet functional` consumes those TSVs directly, which is the same route FRESH's own HUMAnN 3
output will take once available. See [example/functional_demo](../example/functional_demo).
