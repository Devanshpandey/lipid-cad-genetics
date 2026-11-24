#!/bin/bash
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate r_env
cd $SCRATCH/coloc_dist
for L in LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB; do
  for O in CAD MI REVASC STROKE HF CV_DEATH MACE AF PAD; do
    f=out/${L}_vs_${O}_coloc.csv
    [ -s "$f" ] && continue
    Rscript /path/to/cad-genetics/scripts/agent1_genetics/R/run_coloc.R \
      --lipid_trait "$L" --outcome "$O" \
      --sumstat_dir "$SCRATCH/coloc_dist/sumstat_windows" \
      --sentinels   "$SCRATCH/coloc_dist/sentinels_${L}.tsv" \
      --out "$f" >> coloc_runs.log 2>&1
    echo "$L vs $O done ($(ls out/ | wc -l)/54)" >> status.log
  done
done
echo ALL_COLOC_DONE >> status.log
