#!/bin/bash
#------------------------------------------------------------
# Define colocalization sentinels — once per exposure, cached.
# Depends on 07a_coloc_ldpanel.slurm.
# Submit: sbatch --array=0-5 07b_coloc_sentinels.slurm
#------------------------------------------------------------
#SBATCH -J coloc_sent
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 06:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err

set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_R}"
set -u

LIPIDS=(LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB)
TR="${LIPIDS[${SLURM_ARRAY_TASK_ID}]}"

WORKDIR="${SCRATCH}/coloc_rerun"
OUTDIR="${WORKDIR}/sentinels"
mkdir -p "${OUTDIR}"

echo "[$(date '+%H:%M:%S')] Sentinels for ${TR}"

Rscript "${SCRIPTS_DIR}/agent1_genetics/R/make_coloc_sentinels.R" \
  --trait       "${TR}" \
  --sumstat_dir "${AGENT1_OUT}/gwas_summary_stats/merged" \
  --ld_bfile    "${WORKDIR}/coloc_ldref" \
  --plink       "/path/to/software/plink" \
  --clump_r2    0.001 \
  --clump_kb    1000 \
  --min_sep_kb  500 \
  --out         "${OUTDIR}/sentinels_${TR}.tsv"

echo "[$(date '+%H:%M:%S')] ${TR}: $(( $(wc -l < "${OUTDIR}/sentinels_${TR}.tsv") - 1 )) sentinel loci"
echo "[$(date)] Done."
