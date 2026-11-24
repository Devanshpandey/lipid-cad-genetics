#!/bin/bash
#------------------------------------------------------------
# REGENIE Step 1 null model for the TG:HDL ratio (single trait).
# Mirrors 02_regenie_step1.slurm (QT block); reuses existing pruned SNP list.
#------------------------------------------------------------
#SBATCH -J tghdl_step1
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 12:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_REGENIE}"

mkdir -p "${TACC_SCRATCH}/regenie_tmp" "${AGENT1_OUT}/regenie_step1"
PRUNE_SNP_LIST="${AGENT1_OUT}/phenotypes/qc_pruned_snps.txt"
[[ -f "${PRUNE_SNP_LIST}" ]] || { echo "ERROR: pruned SNP list missing: ${PRUNE_SNP_LIST}"; exit 1; }

echo "[$(date)] REGENIE Step 1 — TG:HDL null model"
regenie \
  --step 1 \
  --bed "${UKB_GENO_STEP1}" \
  --extract "${PRUNE_SNP_LIST}" \
  --keep "${AGENT1_OUT}/phenotypes/regenie_keep.txt" \
  --phenoFile "${AGENT1_OUT}/phenotypes/pheno_tghdl.txt" \
  --covarFile "${AGENT1_OUT}/phenotypes/covariates.txt" \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix "${TACC_SCRATCH}/regenie_tmp/step1_tghdl" \
  --threads "${SLURM_NTASKS}" \
  --gz \
  --out "${AGENT1_OUT}/regenie_step1/step1_tghdl"
echo "[$(date)] done -> ${AGENT1_OUT}/regenie_step1/step1_tghdl_pred.list"
