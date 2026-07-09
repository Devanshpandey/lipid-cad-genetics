#!/bin/bash
#------------------------------------------------------------
# PRS-CS: Bayesian continuous shrinkage PRS (auto mode)
# SLURM array: 7 traits × 22 chromosomes = 154 tasks
#   task_id = trait_idx * 22 + (chr - 1)
#   trait_idx: 0=LDL_C, 1=HDL_C, 2=TRIGLY, 3=TOT_CHOL, 4=APOA1, 5=APOB, 6=LPA
#
# Uses UKB EUR LD reference (precomputed hdf5 blocks).
# Input: BETA+SE format (recommended — avoids p<1e-323 truncation).
# Output: per-chr posterior effect sizes → merged by 10b_prs_cs_merge.slurm
#
# Submit: sbatch --array=0-153 10_prs_cs.slurm
#------------------------------------------------------------
#SBATCH -J cad_prscs
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -A OTH26002
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=devanshpandey@utexas.edu

set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/work/07880/devansh/lonestar/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_PY}"
set -u

# PRS-CS is single-threaded MCMC; limit scipy BLAS threads to avoid oversubscription
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

PRSCS_DIR="/work/07880/devansh/lonestar/software/PRScs"
PRSCS_REF="${PRSCS_DIR}/ldblk_ukbb_eur"
BIM_PREFIX="${UKB_GENO_STEP2}"

LIPID_TRAITS=(LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB LPA)
GWAS_N=(354746 325276 355123 355401 323465 353690 282610)

TASK_ID="${SLURM_ARRAY_TASK_ID:-0}"
TRAIT_IDX=$((TASK_ID / 22))
CHR=$((TASK_ID % 22 + 1))
TRAIT="${LIPID_TRAITS[${TRAIT_IDX}]}"
N="${GWAS_N[${TRAIT_IDX}]}"

echo "[$(date)] PRS-CS auto: ${TRAIT} chr${CHR} (task ${TASK_ID}, N=${N})"

OUT_DIR="${AGENT1_OUT}/prs_weights_prscs/${TRAIT}"
mkdir -p "${OUT_DIR}"

# ── Step 1: Prepare sumstat in SNP/A1/A2/BETA/SE format (once per trait) ──
SST_FILE="${OUT_DIR}/${TRAIT}_prscs_sst.txt"
if [[ ! -f "${SST_FILE}" ]]; then
  SUMSTAT="${AGENT1_OUT}/gwas_summary_stats/merged/lipids_${TRAIT}.txt.gz"
  [[ -f "${SUMSTAT}" ]] || { echo "ERROR: GWAS sumstat not found: ${SUMSTAT}"; exit 1; }

  echo "  Preparing sumstat..."
  python3 - <<PYEOF
import pandas as pd
df = pd.read_csv("${SUMSTAT}", sep=r"\s+", na_values=["NA"])
df = df.rename(columns={"ID":"SNP","ALLELE1":"A1","ALLELE0":"A2"})
df = df[["SNP","A1","A2","BETA","SE"]].dropna()
df = df[df["SE"] > 0]
df.to_csv("${SST_FILE}", sep="\t", index=False)
print(f"  {len(df)} SNPs written to ${SST_FILE}")
PYEOF
fi

[[ -f "${SST_FILE}" ]] || { echo "ERROR: sumstat file not created"; exit 1; }

# ── Step 2: Run PRS-CS-auto for this chr ──
OUT_PREFIX="${OUT_DIR}/${TRAIT}_prscs"
EXPECTED="${OUT_PREFIX}_pst_eff_a1_b0.5_phiauto_chr${CHR}.txt"

if [[ -f "${EXPECTED}" ]]; then
  echo "  Output already exists — skipping: ${EXPECTED}"
  echo "  $(wc -l < "${EXPECTED}") SNPs"
else
  echo "  Running PRS-CS-auto..."
  python3 "${PRSCS_DIR}/PRScs.py" \
    --ref_dir="${PRSCS_REF}" \
    --bim_prefix="${BIM_PREFIX}" \
    --sst_file="${SST_FILE}" \
    --n_gwas="${N}" \
    --chrom="${CHR}" \
    --out_dir="${OUT_PREFIX}" \
    --seed=42

  if [[ -f "${EXPECTED}" ]]; then
    echo "  Done: $(wc -l < "${EXPECTED}") SNPs with posterior weights"
  else
    echo "ERROR: PRS-CS output not produced for ${TRAIT} chr${CHR}"
    ls "${OUT_DIR}/"*chr${CHR}* 2>/dev/null
    exit 1
  fi
fi

echo "[$(date)] ${TRAIT} chr${CHR} complete."
