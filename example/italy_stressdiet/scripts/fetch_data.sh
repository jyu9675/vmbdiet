#!/usr/bin/env bash
# Reproducibly fetch the Italian "StressDiet" cohort (Laghi/Foschi et al., Front Cell Infect Microbiol 2025;
# DOI 10.3389/fcimb.2025.1582283; NCBI BioProject PRJNA1188525):
#   - FFQ nutrient table + metabolomics  (journal supplement, CC-BY 4.0)
#   - 16S V3-V4 paired FASTQs            (ENA / SRA)
# Produces: nutrients.tsv, metabolites.tsv, ena_manifest.tsv, reads/*.fastq.gz
set -euo pipefail
cd "$(dirname "$0")"
UA="Mozilla/5.0"; EM="${NCBI_EMAIL:-anonymous@example.org}"
mkdir -p reads suppl

echo "[1/3] journal supplement (nutrients + metabolites)"
curl -sSL -A "$UA" "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12213695/supplementaryFiles" -o suppl.zip
unzip -o suppl.zip -d suppl >/dev/null
python - <<'PY'
import pandas as pd
pd.read_excel("suppl/DataSheet2.xls", engine="xlrd").to_csv("nutrients.tsv", sep="\t", index=False)   # FFQ
pd.read_excel("suppl/DataSheet1.xls", engine="xlrd").to_csv("metabolites.tsv", sep="\t", index=False)  # 1H-NMR
print("  nutrients.tsv, metabolites.tsv written")
PY

echo "[2/3] ENA run manifest (run <-> StressDiet_N <-> fastq URLs)"
curl -sSL -A "$UA" "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1188525&result=read_run&fields=run_accession,sample_alias,read_count,fastq_bytes,fastq_ftp&format=tsv&limit=0" -o ena_runs.tsv
python - <<'PY'
import pandas as pd
t=pd.read_csv("ena_runs.tsv",sep="\t")
t["sample_num"]=t["sample_alias"].str.extract(r"StressDiet_(\d+)").astype(int)
t.sort_values("sample_num").to_csv("ena_manifest.tsv",sep="\t",index=False)
print("  ena_manifest.tsv:",len(t),"runs")
PY

echo "[3/3] download 16S FASTQs (~3 GB)"
awk -F'\t' 'NR>1{n=split($5,a,";"); for(i=1;i<=n;i++) print (a[i]~/^http/?a[i]:"https://"a[i])}' ena_manifest.tsv > urls.txt
sed -i 's/\r$//' urls.txt
xargs -a urls.txt -P 8 -I{} bash -c 'u="$1"; f="reads/$(basename "$u")"; [ -s "$f" ] || curl -sS -L --retry 5 -A "Mozilla/5.0" "$u" -o "$f"' _ {}
echo "done: $(ls reads/*.fastq.gz | wc -l) files"
