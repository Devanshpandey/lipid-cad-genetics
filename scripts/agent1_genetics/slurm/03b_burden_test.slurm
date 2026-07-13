#!/bin/bash
#------------------------------------------------------------
# SLURM array: REGENIE gene-based burden testing (WES)
#
# Data source:
#   ${UKB_WES_PFILE}.pgen/pvar/psam  — all-chr merged pgen, ~200K-470K WES samples
#   (joyce_annotated_vcf/ used by prepare_burden_masks.sh to build annotation files)
#
# Models run per chromosome:
#   A) Quantitative lipid traits — additive
#   B) Quantitative lipid traits — recessive
#   C) Binary CAD/MACE outcomes  — additive (Firth logistic)
#   D) Binary CAD/MACE outcomes  — recessive (Firth logistic)
#
# Two masks per gene:
#   mask1 = pLoF only (LOFTEE HC + high-confidence cryptic splice)
#   mask2 = pLoF + deleterious missense (MiS ≥ 0.7 — as used in Koyama Nat Genet 2026)
#
# AAF bins: 0.001, 0.01 (ultra-rare and rare carrier frequencies)
# Test: ACATO-full (combines SKAT-O + burden tests, optimal for rare variants)
#
# Depends on:
#   Step 02  — REGENIE Step 1 predictions (step1_qt_pred.list, step1_bt_pred.list)
#   prepare_burden_masks.sh — annotation/setlist/mask files must exist
#
# Submit: sbatch --array=1-22 03b_burden_test.slurm
#------------------------------------------------------------
#SBATCH -J cad_burden
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 16:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail
# Resolve repo root: SLURM_SUBMIT_DIR when run via sbatch, hardcoded path as fallback
# (BASH_SOURCE[0] points to /var/spool in SLURM jobs — cannot use dirname approach)
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

CHR="${SLURM_ARRAY_TASK_ID}"

set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_REGENIE}"
set -u

OUT_DIR="${AGENT1_OUT}/burden"
TMP_DIR="${TACC_SCRATCH}/regenie_burden_tmp"
mkdir -p "${OUT_DIR}" "${TMP_DIR}"

echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME:-shell}/${SLURM_ARRAY_JOB_ID:-0}_${CHR}] REGENIE burden test — Chr ${CHR}"
echo "Threads:  ${SLURM_NTASKS}"
echo "WES pfile: ${UKB_WES_PFILE}.{pgen,pvar,psam}"

# ---- Preflight checks ----
for ext in pgen pvar psam; do
  [[ -f "${UKB_WES_PFILE}.${ext}" ]] || {
    echo "ERROR: WES pfile not found: ${UKB_WES_PFILE}.${ext}"
    echo "Check UKB_WES_PFILE in config/tacc_paths.sh"
    exit 1
  }
done
for f in "${UKB_BURDEN_ANNO_FILE}" "${UKB_BURDEN_SETLIST}" "${UKB_BURDEN_MASKS}"; do
  [[ -f "${f}" ]] || {
    echo "ERROR: Burden annotation file missing: ${f}"
    echo "Run: bash scripts/agent1_genetics/slurm/prepare_burden_masks.sh"
    exit 1
  }
done

STEP1_QT_PRED="${AGENT1_OUT}/regenie_step1/step1_qt_pred.list"
STEP1_BT_PRED="${AGENT1_OUT}/regenie_step1/step1_bt_pred.list"
KEEP_FILE="${AGENT1_OUT}/phenotypes/regenie_keep.txt"
PHENO_QT="${AGENT1_OUT}/phenotypes/pheno_quantitative.txt"
PHENO_BT="${AGENT1_OUT}/phenotypes/pheno_binary.txt"
COVAR="${AGENT1_OUT}/phenotypes/covariates.txt"

for f in "${STEP1_QT_PRED}" "${STEP1_BT_PRED}" "${KEEP_FILE}" "${PHENO_QT}" "${PHENO_BT}" "${COVAR}"; do
  [[ -f "${f}" ]] || { echo "ERROR: Required input missing: ${f}"; exit 1; }
done

# ---- Common REGENIE flags for all 4 runs ----
# Using --pgen (PLINK2 pgen format) + --chr to restrict to one chromosome.
# --keep intersects WES samples with white-British GWAS samples.
# --minMAC 1 keeps singletons (needed for ultra-rare burden — REGENIE handles AC in masks).
# NOTE: --vc-tests and --test are mutually exclusive in REGENIE v4.
#   Additive runs use --vc-tests acato-full (omnibus: SKAT+burden+ACAT).
#   Recessive runs drop --vc-tests and use --test recessive (collapsing burden).
COMMON_FLAGS=(
  --pgen    "${UKB_WES_PFILE}"
  --keep    "${KEEP_FILE}"
  --covarFile "${COVAR}"
  --anno-file "${UKB_BURDEN_ANNO_FILE}"
  --set-list  "${UKB_BURDEN_SETLIST}"
  --mask-def  "${UKB_BURDEN_MASKS}"
  --aaf-bins  0.001,0.01
  --build-mask max
  --chr       "${CHR}"
  --minMAC    1
  --threads   "${SLURM_NTASKS}"
  --gz
)
ADDITIVE_FLAGS=( "${COMMON_FLAGS[@]}" --vc-tests acato-full )
RECESSIVE_FLAGS=( "${COMMON_FLAGS[@]}" --test recessive )

# ---- A: Quantitative traits — additive ----
OUT_A="${OUT_DIR}/qt_additive_chr${CHR}"
if ls "${OUT_A}"_*.regenie.gz &>/dev/null 2>&1; then
  echo "[$(date)] A (QT additive chr${CHR}) already done — skipping"
else
  echo "[$(date)] A: Quantitative additive..."
  regenie \
    --step 2 \
    "${ADDITIVE_FLAGS[@]}" \
    --phenoFile "${PHENO_QT}" \
    --pred      "${STEP1_QT_PRED}" \
    --apply-rint \
    --lowmem --lowmem-prefix "${TMP_DIR}/burden_qt_add_chr${CHR}" \
    --out "${OUT_A}"
  echo "[$(date)] A complete: ${OUT_A}*.regenie.gz"
fi

# ---- B: Quantitative traits — recessive ----
# Koyama Nat Genet 2026: recessive effects for ANGPTL4, TM6SF2, BTNL9, ABHD15
# are 2-3× larger than additive — these are missed by additive-only analysis.
OUT_B="${OUT_DIR}/qt_recessive_chr${CHR}"
if ls "${OUT_B}"_*.regenie.gz &>/dev/null 2>&1; then
  echo "[$(date)] B (QT recessive chr${CHR}) already done — skipping"
else
  echo "[$(date)] B: Quantitative recessive..."
  regenie \
    --step 2 \
    "${RECESSIVE_FLAGS[@]}" \
    --phenoFile "${PHENO_QT}" \
    --pred      "${STEP1_QT_PRED}" \
    --apply-rint \
    --lowmem --lowmem-prefix "${TMP_DIR}/burden_qt_rec_chr${CHR}" \
    --out "${OUT_B}"
  echo "[$(date)] B complete: ${OUT_B}*.regenie.gz"
fi

# ---- C: Binary traits — additive (Firth logistic) ----
OUT_C="${OUT_DIR}/bt_additive_chr${CHR}"
if ls "${OUT_C}"_*.regenie.gz &>/dev/null 2>&1; then
  echo "[$(date)] C (BT additive chr${CHR}) already done — skipping"
else
  echo "[$(date)] C: Binary additive..."
  regenie \
    --step 2 \
    "${ADDITIVE_FLAGS[@]}" \
    --phenoFile "${PHENO_BT}" \
    --pred      "${STEP1_BT_PRED}" \
    --bt --firth --approx \
    --out "${OUT_C}"
  echo "[$(date)] C complete: ${OUT_C}*.regenie.gz"
fi

# ---- D: Binary traits — recessive (Firth logistic) ----
OUT_D="${OUT_DIR}/bt_recessive_chr${CHR}"
if ls "${OUT_D}"_*.regenie.gz &>/dev/null 2>&1; then
  echo "[$(date)] D (BT recessive chr${CHR}) already done — skipping"
else
  echo "[$(date)] D: Binary recessive..."
  regenie \
    --step 2 \
    "${RECESSIVE_FLAGS[@]}" \
    --phenoFile "${PHENO_BT}" \
    --pred      "${STEP1_BT_PRED}" \
    --bt --firth --approx \
    --out "${OUT_D}"
  echo "[$(date)] D complete: ${OUT_D}*.regenie.gz"
fi

echo ""
echo "[$(date)] Chr ${CHR} — all 4 burden runs complete."
echo "Outputs in ${OUT_DIR}/*chr${CHR}*.regenie.gz"
