#!/bin/bash
#SBATCH -J statbp
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/statbp_%j.out
#SBATCH -e /path/to/cad-genetics/statbp_%j.err
# Reviewer round 4, item 3: strengthen the prevalent-CAD exclusion with baseline
# self-report and OPCS-4 revascularization, then report whether the score-by-statin
# interactions still hold. Block 3 of 08l was invoked with --out alone and could
# never have run; the score/weight/raw/covar arguments below are the ones 08e uses.
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
W=$SCRATCH/statint
for f in $W/dose_LPA.raw $W/dose_LDL.raw $W/lpa_score.txt $W/ldl_score.txt; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f -- rerun 08e's score build first"; exit 1; }
done
Rscript $_REPO/scripts/agent1_genetics/R/statin_interaction_incident_broadprev.R \
  --scores  $W/dose_LPA.raw,$W/dose_LDL.raw \
  --weights $W/lpa_score.txt,$W/ldl_score.txt \
  --raw_dir /path/to/ukb_phenotypes \
  --covar   $_REPO/results/agent1_genetics/phenotypes/covariates.txt \
  --out     $_REPO/exports/agent1/revision/statin_incident_broadprev.csv || echo "BROADPREV FAILED"
echo ALL_DONE
