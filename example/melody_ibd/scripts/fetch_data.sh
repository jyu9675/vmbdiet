#!/usr/bin/env bash
# Reproducibly fetch the MELODY vaginal cohort (PLOS ONE 2024, "Association of vaginal cytokines
# and dietary intake with IBD status and vaginal microbiota in pregnant individuals";
# DOI 10.1371/journal.pone.0335178; BioProject PRJNA915128):
#   - vaginal 16S V3-V4 FASTQs (ENA; the project also holds unrelated stool WGS runs, excluded here)
#   - per-sample diet (HEI-2015 + nutrients) + covariates, parsed from the BioSample metadata
# Produces: manifest.tsv, run_alias.tsv, diet.tsv, reads/*.fastq.gz
set -euo pipefail
cd "$(dirname "$0")"; UA="Mozilla/5.0"; mkdir -p reads samplexml

echo "[1/3] ENA manifest -> vaginal AMPLICON subset only"
curl -sSL -A "$UA" "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA915128&result=read_run&fields=run_accession,sample_accession,sample_alias,scientific_name,library_strategy,fastq_ftp&format=tsv&limit=0" -o ena_all.tsv
python - <<'PY'
import pandas as pd
t=pd.read_csv("ena_all.tsv",sep="\t")
v=t[(t.scientific_name=="vaginal metagenome")&(t.library_strategy=="AMPLICON")].copy()
v.to_csv("manifest.tsv",sep="\t",index=False)
v[["run_accession","sample_alias"]].to_csv("run_alias.tsv",sep="\t",index=False,lineterminator="\n")
urls=[("https://"+u) if not u.startswith("http") else u
      for b in v.fastq_ftp.dropna() for u in str(b).split(";") if u]
open("urls.txt","w",newline="\n").write("\n".join(urls)+"\n")
print("  vaginal amplicon runs:",len(v))
PY

echo "[2/3] diet + covariates from BioSample metadata"
cut -f2 manifest.tsv | tail -n +2 | sort -u > samples.txt
cat samples.txt | xargs -P 8 -I{} bash -c 'curl -sS -L -A "Mozilla/5.0" "https://www.ebi.ac.uk/ena/browser/api/xml/{}" -o "samplexml/{}.xml"'
python - <<'PY'
import re,glob,pandas as pd
def A(fp):
    x=open(fp,encoding="utf-8",errors="ignore").read()
    d={t.strip():v.strip() for t,v in re.findall(r"<TAG>(.*?)</TAG>\s*<VALUE>(.*?)</VALUE>",x,re.S)}
    m=re.search(r'accession="(SAM[NED][A-Z]?\d+)"',x); d["_s"]=m.group(1) if m else fp
    return d
def g(d,*k):
    for c in k:
        if c in d: return d[c]
    return None
rows=[]
for fp in glob.glob("samplexml/*.xml"):
    d=A(fp)
    rows.append(dict(sample=d["_s"], sample_alias=g(d,"id_merged","id_merged_simple"),
      subject=g(d,"study_id_mom","study_id_mom1"), ibd=g(d,"cd_uc_hc_diagnosis_x","ibd_group_x"),
      age=g(d,"age"), bmi=g(d,"bmi"), race=g(d,"race"),
      HEI2015=g(d,"hei_2015_total_score"), AHEI=g(d,"ahei_10_score_0_110_","ahei_10_score"),
      animal_protein=g(d,"animal_protein"), vegetable_protein=g(d,"vegetable_protein"),
      total_fat=g(d,"total_fat"), saturated_fat=g(d,"total_saturated_fatty_acids_sfa_"),
      fiber=g(d,"total_dietary_fiber"), carbohydrate=g(d,"total_carbohydrate"),
      alcohol=g(d,"alcohol"), energy_kcal=g(d,"energy")))
m=pd.DataFrame(rows)
for c in ["age","bmi","HEI2015","AHEI","animal_protein","vegetable_protein","total_fat",
          "saturated_fat","fiber","carbohydrate","alcohol","energy_kcal"]:
    m[c]=pd.to_numeric(m[c],errors="coerce")
m.sort_values("sample_alias").to_csv("diet.tsv",sep="\t",index=False,lineterminator="\n")
print("  diet.tsv:",len(m),"samples")
PY

echo "[3/3] download vaginal FASTQs (~1.3 GB)"
sed -i 's/\r$//' urls.txt
xargs -a urls.txt -P 8 -I{} bash -c 'u="$1"; f="reads/$(basename "$u")"; [ -s "$f" ] || curl -sS -L --retry 5 -A "Mozilla/5.0" "$u" -o "$f"' _ {}
echo "done: $(ls reads/*.fastq.gz | wc -l) files"
