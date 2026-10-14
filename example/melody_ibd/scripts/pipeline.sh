#!/usr/bin/env bash
# 16S V3-V4 processing: MELODY (PRJNA915128, vaginal subset) -> per-sample ASV/genus+species table.
# Primers already removed in the deposited reads, so NO cutadapt step (unlike the Italian cohort):
# truncate -> vsearch merge -> EE filter -> UNOISE3 ASVs -> chimera -> otutab -> SINTAX(+species).
set -euo pipefail
# Portable: BASE is this example dir. Provide tools via env vars:
#   VSEARCH  path to vsearch executable (default: `vsearch` on PATH)
#   RDP_REF  path to the SINTAX RDP reference (default: ref/rdp_16s_v16.fa.gz under BASE)
BASE="$(cd "$(dirname "$0")/.." && pwd)"
cd "$BASE"
VS="${VSEARCH:-vsearch}"
REF="${RDP_REF:-$BASE/ref/rdp_16s_v16.fa.gz}"
mkdir -p work/filt work/tmp logs

echo "[$(date +%T)] per-sample truncate/merge/filter (no primer trim)"
MAP="$BASE/data/run_alias.tsv"; [ -f "$MAP" ] || MAP="$BASE/run_alias.tsv"
tail -n +2 "$MAP" | while IFS=$'\t' read -r srr alias; do
  r1="reads/${srr}_1.fastq.gz"; r2="reads/${srr}_2.fastq.gz"; out="work/filt/${alias}.fa"
  [ -s "$out" ] && { echo "  skip $alias"; continue; }
  [ -s "$r1" ] && [ -s "$r2" ] || { echo "  MISSING $srr ($alias)"; continue; }
  t=work/tmp
  "$VS" --fastx_filter "$r1" --fastq_trunclen_keep 270 --fastqout $t/t1.fq --quiet
  "$VS" --fastx_filter "$r2" --fastq_trunclen_keep 210 --fastqout $t/t2.fq --quiet
  "$VS" --fastq_mergepairs $t/t1.fq --reverse $t/t2.fq --fastqout $t/m.fq \
     --fastq_maxdiffs 15 --fastq_minovlen 20 --quiet
  "$VS" --fastx_filter $t/m.fq --fastq_maxee 1.0 --fastq_minlen 380 --fastq_maxlen 460 \
     --relabel "${alias}." --fastaout "$out" --quiet
  echo "  $alias  $(grep -c '^>' "$out") seqs"
done

echo "[$(date +%T)] pool + dereplicate"
cat work/filt/*.fa > work/all.fa
grep -c '^>' work/all.fa | xargs echo "  total reads:"
"$VS" --derep_fulllength work/all.fa --sizeout --minuniquesize 2 --output work/uniques.fa --quiet
echo "  uniques: $(grep -c '^>' work/uniques.fa)"
echo "[$(date +%T)] UNOISE3 + chimera"
"$VS" --cluster_unoise work/uniques.fa --minsize 8 --sizein --sizeout --centroids work/zc.fa --quiet
"$VS" --uchime3_denovo work/zc.fa --sizein --nonchimeras work/zotus.fa --relabel ASV --quiet
echo "  ASVs: $(grep -c '^>' work/zotus.fa)"
echo "[$(date +%T)] OTU table"
"$VS" --usearch_global work/all.fa --db work/zotus.fa --id 0.97 --strand plus --otutabout work/otutab.txt --quiet
echo "[$(date +%T)] SINTAX"
"$VS" --sintax work/zotus.fa --db "$REF" --tabbedout work/sintax.txt --strand both --sintax_cutoff 0.8 --quiet
echo "[$(date +%T)] DONE pipeline (ASVs=$(grep -c '^>' work/zotus.fa), samples=$(($(head -1 work/otutab.txt|awk '{print NF}')-1)))"
