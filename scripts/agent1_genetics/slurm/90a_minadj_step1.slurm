#!/bin/bash
#SBATCH -J minadj_s1
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
_REPO="/path/to/cad-genetics"
source "${_REPO}/config/tacc_paths.sh"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_REGENIE}"
mkdir -p "${TACC_SCRATCH}/regenie_tmp" "${AGENT1_OUT}/regenie_step1"
echo "[$(date)] REGENIE Step 1 — NMR LDL particle traits"
regenie \
  --step 1 \
  --bed "${UKB_GENO_STEP1}" \
  --extract "${AGENT1_OUT}/phenotypes/qc_pruned_snps.txt" \
  --keep "${AGENT1_OUT}/phenotypes/regenie_keep.txt" \
  --phenoFile "${AGENT1_OUT}/phenotypes/pheno_lipids.txt" \
  --covarFile "${AGENT1_OUT}/phenotypes/covariates_min.txt" \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix "${TACC_SCRATCH}/regenie_tmp/step1_minadj" \
  --threads "${SLURM_NTASKS}" \
  --gz \
  --out "${AGENT1_OUT}/regenie_step1/step1_minadj"
echo "[$(date)] NMR Step 1 complete"
