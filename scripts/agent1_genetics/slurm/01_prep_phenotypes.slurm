#!/bin/bash
#------------------------------------------------------------
# SLURM: Phenotype preparation (runs on login or compute node)
#------------------------------------------------------------
#SBATCH -J cad_pheno_prep
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 02:00:00
#SBATCH -A OTH26002
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=devanshpandey@utexas.edu

set -euo pipefail
# Resolve repo root: SLURM_SUBMIT_DIR when run via sbatch, hardcoded path as fallback
# (BASH_SOURCE[0] points to /var/spool in SLURM jobs — cannot use dirname approach)
_REPO="${SLURM_SUBMIT_DIR:-/work/07880/devansh/lonestar/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_R}"

mkdir -p "${AGENT1_OUT}/phenotypes"

echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME:-shell}/${SLURM_JOB_ID:-0}] Starting phenotype preparation"

Rscript "${SCRIPTS_DIR}/agent1_genetics/R/prep_phenotypes.R" \
  --pheno_dir "${UKB_PHENO_DIR}" \
  --icd_file  "${UKB_ICD_FILE}" \
  --sample    "${UKB_SAMPLE_EID}" \
  --out       "${AGENT1_OUT}/phenotypes"

echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME:-shell}/${SLURM_JOB_ID:-0}] Phenotype preparation complete"
echo "Outputs in: ${AGENT1_OUT}/phenotypes/"
