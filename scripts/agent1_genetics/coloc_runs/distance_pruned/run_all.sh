#!/bin/bash
source /path/to/anaconda3/etc/profile.d/conda.sh
cd $SCRATCH/coloc_dist
conda activate pyt_env
python3 extract_windows.py > extract.log 2>&1 || { echo EXTRACT_FAILED >> status.log; exit 1; }
echo EXTRACT_OK >> status.log
conda activate r_env
mkdir -p out
LIP="LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB"
OUTC="CAD MI REVASC STROKE HF CV_DEATH MACE AF PAD"
n=0
for L in $LIP; do
  for O in $OUTC; do
    f=out/${L}_vs_${O}_coloc.csv
    [ -s "$f" ] && continue
    Rscript /path/to/cad-genetics/scripts/agent1_genetics/R/run_coloc.R \
      --lipid_trait "$L" --outcome "$O" \
      --sumstat_dir "$SCRATCH/coloc_dist/sumstat_windows" \
      --sentinels   "$SCRATCH/coloc_dist/sentinels_${L}.tsv" \
      --out "$f" >> coloc_runs.log 2>&1
    n=$((n+1)); echo "$L vs $O done ($n/54)" >> status.log
  done
done
python3 - <<'PY' >> status.log 2>&1
import glob,csv,os
rows=[]
for f in glob.glob(os.path.expandvars("$SCRATCH/coloc_dist/out/*_coloc.csv")):
    rows += list(csv.DictReader(open(f)))
if rows:
    w=csv.DictWriter(open(os.path.expandvars("$SCRATCH/coloc_dist/coloc_results_v2.csv"),"w"),fieldnames=rows[0].keys())
    w.writeheader(); w.writerows(rows)
    sig=[r for r in rows if float(r["PP.H4"])>0.8]
    print(f"COMBINED {len(rows)} tests, {len(sig)} at PP.H4>0.8")
PY
echo ALL_COLOC_DONE >> status.log
