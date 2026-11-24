#!/bin/bash
source /path/to/anaconda3/etc/profile.d/conda.sh; conda activate r_env
cd $SCRATCH/coloc_dist
/usr/bin/time -f "patched runtime: %E" Rscript /path/to/cad-genetics/scripts/agent1_genetics/R/run_coloc.R \
  --lipid_trait LDL_C --outcome CAD \
  --sumstat_dir "$SCRATCH/coloc_dist/sumstat_windows" \
  --sentinels "$SCRATCH/coloc_dist/sentinels_LDL_C.tsv" \
  --out validate_LDL_C_vs_CAD.csv > validate.log 2>&1
echo VALIDATE_DONE >> validate.log
