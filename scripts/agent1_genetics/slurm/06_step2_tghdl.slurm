#!/bin/bash
#------------------------------------------------------------
# REGENIE gene-burden test for TG:HDL ratio (additive, ACATO-full).
# Mirrors 03b_burden_test.slurm block A. Submit: sbatch --array=1-22 06_step2_tghdl.slurm
#------------------------------------------------------------
#SBATCH -J tghdl_burden
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 8:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"
CHR="${SLURM_ARRAY_TASK_ID}"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_REGENIE}"; set -u

OUT_DIR="${AGENT1_OUT}/burden"; TMP_DIR="${TACC_SCRATCH}/regenie_burden_tmp"
mkdir -p "${OUT_DIR}" "${TMP_DIR}"
PRED="${AGENT1_OUT}/regenie_step1/step1_tghdl_pred.list"
PHENO="${AGENT1_OUT}/phenotypes/pheno_tghdl.txt"
for f in "${PRED}" "${PHENO}" "${UKB_BURDEN_ANNO_FILE}" "${UKB_BURDEN_SETLIST}" "${UKB_BURDEN_MASKS}"; do
  [[ -f "${f}" ]] || { echo "ERROR: missing ${f}"; exit 1; }
done

OUT_A="${OUT_DIR}/tghdl_additive_chr${CHR}"
if ls "${OUT_A}"_*.regenie.gz &>/dev/null 2>&1; then
  echo "[$(date)] chr${CHR} already done — skipping"
else
  echo "[$(date)] TG:HDL burden additive — chr ${CHR}"
  regenie \
    --step 2 \
    --pgen "${UKB_WES_PFILE}" \
    --keep "${AGENT1_OUT}/phenotypes/regenie_keep.txt" \
    --covarFile "${AGENT1_OUT}/phenotypes/covariates.txt" \
    --anno-file "${UKB_BURDEN_ANNO_FILE}" \
    --set-list  "${UKB_BURDEN_SETLIST}" \
    --mask-def  "${UKB_BURDEN_MASKS}" \
    --aaf-bins  0.001,0.01 \
    --build-mask max \
    --chr "${CHR}" \
    --minMAC 1 \
    --threads "${SLURM_NTASKS}" \
    --gz \
    --vc-tests acato-full \
    --phenoFile "${PHENO}" \
    --pred "${PRED}" \
    --apply-rint \
    --lowmem --lowmem-prefix "${TMP_DIR}/burden_tghdl_chr${CHR}" \
    --out "${OUT_A}"
  echo "[$(date)] chr${CHR} complete: ${OUT_A}*.regenie.gz"
fi
