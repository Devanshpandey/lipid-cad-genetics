#!/bin/bash
source /path/to/anaconda3/etc/profile.d/conda.sh; conda activate r_env
_REPO=/path/to/cad-genetics; W=$SCRATCH/statint
Rscript $_REPO/scripts/agent1_genetics/R/statin_interaction_robustness.R \
  --scores  $W/dose_LPA.raw,$W/dose_LDL.raw \
  --weights $W/lpa_score.txt,$W/ldl_score.txt \
  --raw_dir /path/to/ukb_phenotypes \
  --covar   /path/to/cad-genetics/results/agent1_genetics/phenotypes/covariates.txt \
  --out     $W/statin_robustness_fixed.csv > $W/rerun.log 2>&1
echo RERUN_DONE >> $W/rerun.log
