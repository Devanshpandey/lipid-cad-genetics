#!/bin/bash
#------------------------------------------------------------
# Cross-trait coloc on the LD-CLUMPED sentinel definition.
#
# Companion to $SCRATCH/coloc_dist (distance-pruned sentinels). The ONLY
# difference is --sentinels: same run_coloc.R, same +/-500 kb windows, same
# priors, and deliberately the SAME --sumstat_dir, because every clumped
# sentinel was verified to be present in those extracted windows (1,592/1,592).
# Sharing the inputs byte-for-byte means any difference in the results is
# attributable to the sentinel definition alone and to nothing else.
#
# Array index = one lipid; each task loops its 9 outcomes.
# Submit: sbatch --array=0-5 run_clump_coloc.slurm
#------------------------------------------------------------
#SBATCH -J coloc_clump
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o coloc_clump_%A_%a.out
#SBATCH -e coloc_clump_%A_%a.err
#SBATCH --mail-type=fail
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail
source /path/to/anaconda3/etc/profile.d/conda.sh
conda activate r_env
cd $SCRATCH/coloc_clump

LIPIDS=(LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB)
L=${LIPIDS[$SLURM_ARRAY_TASK_ID]}
WIN=$SCRATCH/coloc_dist/sumstat_windows
SENT=$SCRATCH/coloc_clump/sentinels_${L}.tsv

[ -s "$SENT" ] || { echo "FATAL: missing sentinels $SENT"; exit 1; }
[ -d "$WIN" ]  || { echo "FATAL: missing window dir $WIN"; exit 1; }

for O in CAD MI REVASC STROKE HF CV_DEATH MACE AF PAD; do
  f=out/${L}_vs_${O}_coloc.csv
  [ -s "$f" ] && { echo "skip $L vs $O (exists)"; continue; }
  echo "=== $L vs $O ==="
  Rscript /path/to/cad-genetics/scripts/agent1_genetics/R/run_coloc.R \
    --lipid_trait "$L" --outcome "$O" \
    --sumstat_dir "$WIN" \
    --sentinels   "$SENT" \
    --out "$f"
done
echo "TASK_DONE $L"
