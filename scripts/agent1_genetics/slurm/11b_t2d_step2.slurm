#!/bin/bash
#------------------------------------------------------------
# REGENIE Step 2 for T2D: (A) chr11 gene-burden — the decisive PDE3B LoF test
# — and (B) chr11 common variants, for cis-MR and coloc of the PDE3B TG signal
# against T2D. Depends on 11a_t2d_step1.slurm.
#
# Burden masks are read from the corral archive (they were cleared from /work).
# Submit: sbatch --dependency=afterok:<step1 jobid> 11b_t2d_step2.slurm
#------------------------------------------------------------
#SBATCH -J t2d_step2
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 12:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_REGENIE}"; set -u

W="${SCRATCH}/t2d"
ARCH=/path/to/archive/results/agent1_genetics/burden_masks
OUT="${W}/step2"; mkdir -p "${OUT}"
PRED="${W}/step1/step1_t2d_pred.list"
[ -s "${PRED}" ] || { echo "FATAL: step1 predictions missing: ${PRED}"; exit 1; }

# ---- A: gene burden, chr11 (PDE3B) ----
echo "[$(date)] A: chr11 gene burden for T2D"
regenie \
  --step 2 \
  --pgen "${UKB_WES_PFILE}" \
  --keep "${AGENT1_OUT}/phenotypes/regenie_keep.txt" \
  --phenoFile "${W}/pheno_t2d.txt" \
  --covarFile "${W}/covariates_t2d.txt" \
  --pred "${PRED}" \
  --anno-file "${ARCH}/ukb_wes_anno.txt" \
  --set-list  "${ARCH}/ukb_wes_setlist.txt" \
  --mask-def  "${ARCH}/ukb_wes_masks.txt" \
  --aaf-bins 0.001,0.01 \
  --build-mask max \
  --chr 11 \
  --minMAC 1 \
  --bt --firth --approx \
  --vc-tests acato-full \
  --threads "${SLURM_NTASKS}" \
  --gz \
  --out "${OUT}/t2d_burden_chr11"

echo "[$(date)] PDE3B burden result:"
zcat "${OUT}"/t2d_burden_chr11_*.regenie.gz 2>/dev/null | grep PDE3B | head -10

# ---- B: chr11 common variants, for cis-MR / coloc at PDE3B ----
echo "[$(date)] B: chr11 common-variant scan for T2D"
regenie \
  --step 2 \
  --bed "${UKB_GENO_STEP2}" \
  --keep "${AGENT1_OUT}/phenotypes/regenie_keep.txt" \
  --phenoFile "${W}/pheno_t2d.txt" \
  --covarFile "${W}/covariates_t2d.txt" \
  --pred "${PRED}" \
  --chr 11 \
  --bsize 400 \
  --minMAC 20 \
  --bt --firth --approx \
  --threads "${SLURM_NTASKS}" \
  --gz \
  --out "${OUT}/t2d_common_chr11"

echo "[$(date)] Done."
ls -la "${OUT}"
