#!/bin/bash
#------------------------------------------------------------
# REGENIE Step 1 for type-2 diabetes — PDE3B target-safety analysis.
#
# WHY. PDE3B loss-of-function lowers triglycerides by 0.42 SD and raises HDL-C
# by 0.35 SD in UKB exomes (~1,020 carriers), the largest lipid effect of any
# drug target in our set. PDE3B mediates insulin's antilipolytic action in
# adipocytes, so the obvious safety question is whether that TG benefit costs
# glycaemic control. The existing HbA1c burden result is null (-0.028 SD,
# 95% CI -0.079 to +0.024) but HbA1c is treated down in diagnosed diabetics,
# which biases toward exactly that null. A binary T2D endpoint avoids it.
#
# Covariates deliberately EXCLUDE `diabetes` (it defines the outcome) and
# `bmi` (a mediator of the PDE3B->T2D path). A BMI-adjusted covariate file is
# written alongside for sensitivity only.
#
# Submit: sbatch 11a_t2d_step1.slurm
#------------------------------------------------------------
#SBATCH -J t2d_step1
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 24:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_REGENIE}"; set -u

W="${SCRATCH}/t2d"
mkdir -p "${W}/step1" "${SCRATCH}/regenie_tmp"

echo "[$(date)] T2D step 1 — $(wc -l < ${W}/pheno_t2d.txt) samples with a T2D call"

regenie \
  --step 1 \
  --bed "${UKB_GENO_STEP1}" \
  --extract "${AGENT1_OUT}/phenotypes/qc_pruned_snps.txt" \
  --keep "${AGENT1_OUT}/phenotypes/regenie_keep.txt" \
  --phenoFile "${W}/pheno_t2d.txt" \
  --covarFile "${W}/covariates_t2d.txt" \
  --bsize 1000 \
  --bt \
  --lowmem \
  --lowmem-prefix "${SCRATCH}/regenie_tmp/step1_t2d" \
  --threads "${SLURM_NTASKS}" \
  --gz \
  --out "${W}/step1/step1_t2d"

echo "[$(date)] Step 1 complete:"
ls -la "${W}/step1/"
