#!/usr/bin/env bash
# 16S V3-V4 processing: PRJNA1188525 (StressDiet, Italy) -> per-sample ASV/genus table
# cutadapt (primer trim) -> truncate -> vsearch merge -> EE filter -> UNOISE3 ASVs -> chimera -> otutab -> SINTAX
set -euo pipefail
BASE="/c/Users/dotne/AppData/Local/Temp/claude/C--Users-dotne/f6d0e8d6-d4bf-47ac-9c34-9bbe1ab7c80c/scratchpad/vmbdiet_italy"
cd "$BASE"
VS="$BASE/tools/vsearch-2.32.0-win-x86_64/bin/vsearch.exe"
REF="$BASE/ref/rdp_16s_v16.fa.gz"
FWD="CCTACGGGNGGCWGCAG"; REV="GACTACHVGGGTATCTAATCC"
mkdir -p work/filt work/tmp logs
MAN="$BASE/ena_manifest.tsv"

echo "[$(date +%T)] per-sample trim/merge/filter"
# manifest cols: run_accession sample_alias read_count fastq_bytes fastq_ftp sample_num
tail -n +2 "$MAN" | while IFS=$'\t' read -r srr alias rc fb ftp snum; do
  r1="reads/${srr}_1.fastq.gz"; r2="reads/${srr}_2.fastq.gz"
  out="work/filt/${alias}.fa"
  [ -s "$out" ] && { echo "  skip $alias"; continue; }
  [ -s "$r1" ] && [ -s "$r2" ] || { echo "  MISSING reads $srr ($alias)"; continue; }
  t=work/tmp
  python -m cutadapt -g "$FWD" -G "$REV" --discard-untrimmed -j 0 \
     -o $t/c1.fq.gz -p $t/c2.fq.gz "$r1" "$r2" > logs/cutadapt_${alias}.log 2>&1
  "$VS" --fastx_filter $t/c1.fq.gz --fastq_trunclen_keep 270 --fastqout $t/t1.fq --quiet
  "$VS" --fastx_filter $t/c2.fq.gz --fastq_trunclen_keep 210 --fastqout $t/t2.fq --quiet
  "$VS" --fastq_mergepairs $t/t1.fq --reverse $t/t2.fq --fastqout $t/m.fq \
     --fastq_maxdiffs 15 --fastq_minovlen 20 --quiet
  "$VS" --fastx_filter $t/m.fq --fastq_maxee 1.0 --fastq_minlen 380 --fastq_maxlen 450 \
     --relabel "${alias}." --fastaout "$out" --quiet
  echo "  $alias  $(grep -c '^>' "$out") seqs"
done

echo "[$(date +%T)] pooling + dereplicate"
cat work/filt/*.fa > work/all.fa
grep -c '^>' work/all.fa | xargs echo "  total merged/filtered reads:"
"$VS" --derep_fulllength work/all.fa --sizeout --minuniquesize 2 --output work/uniques.fa --quiet
echo "  uniques: $(grep -c '^>' work/uniques.fa)"

echo "[$(date +%T)] UNOISE3 denoise -> ASVs + chimera removal"
"$VS" --cluster_unoise work/uniques.fa --minsize 8 --sizein --sizeout \
   --centroids work/zotus_chim.fa --quiet
"$VS" --uchime3_denovo work/zotus_chim.fa --sizein --nonchimeras work/zotus.fa \
   --relabel ASV --quiet
echo "  ASVs: $(grep -c '^>' work/zotus.fa)"

echo "[$(date +%T)] OTU table (map reads to ASVs)"
"$VS" --usearch_global work/all.fa --db work/zotus.fa --id 0.97 --strand plus \
   --otutabout work/otutab.txt --quiet
echo "  otutab lines: $(wc -l < work/otutab.txt)"

echo "[$(date +%T)] SINTAX taxonomy (RDP, genus)"
"$VS" --sintax work/zotus.fa --db "$REF" --tabbedout work/sintax.txt \
   --strand both --sintax_cutoff 0.8 --quiet
echo "  sintax lines: $(wc -l < work/sintax.txt)"
echo "[$(date +%T)] DONE pipeline"
