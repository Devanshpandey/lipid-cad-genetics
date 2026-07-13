#!/bin/bash
#------------------------------------------------------------
# PRS: Build PRS weights for each lipid trait (C+T approach)
# SLURM array: one trait per task (0=LDL_C, 1=HDL_C, ..., 6=LPA)
#
# Clumping: PLINK1.9 --clump at GWS (p < 5e-8, r2 < 0.001, 10Mb)
# Scoring:  PLINK1.9 --score for all UKB individuals
#
# Submit: sbatch --array=0-6 10_prs.slurm
#------------------------------------------------------------
#SBATCH -J cad_prs
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 08:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

PLINK1="/path/to/software/plink"
[[ -f "${PLINK1}" ]] || { echo "ERROR: plink1.9 not found at ${PLINK1}" >&2; exit 1; }

set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_PY}"
set -u

LIPID_TRAITS=(LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB LPA)
TRAIT="${LIPID_TRAITS[${SLURM_ARRAY_TASK_ID:-0}]}"
N_THREADS="${SLURM_NTASKS:-16}"

echo "[$(date)] === PRS for ${TRAIT} (task ${SLURM_ARRAY_TASK_ID:-0}) ==="

OUT_DIR="${AGENT1_OUT}/prs_weights/${TRAIT}"
mkdir -p "${OUT_DIR}"

SUMSTAT="${AGENT1_OUT}/gwas_summary_stats/merged/lipids_${TRAIT}.txt.gz"
[[ -f "${SUMSTAT}" ]] || { echo "ERROR: sumstat not found: ${SUMSTAT}"; exit 1; }

PRSICE_INPUT="${OUT_DIR}/${TRAIT}_prsice_input.txt"

# ── Step 1: Prepare sumstat ──
if [[ ! -f "${PRSICE_INPUT}" ]]; then
  echo "  Preparing sumstat..."
  python3 - <<PYEOF
import pandas as pd
df = pd.read_csv("${SUMSTAT}", sep=r"\s+", na_values=["NA"])
df["P"] = 10.0 ** (-df["LOG10P"])
df = df.rename(columns={"ID":"SNP","ALLELE1":"A1","ALLELE0":"A2","A1FREQ":"MAF"})
df = df[["SNP","A1","A2","BETA","SE","P","MAF"]].dropna(subset=["BETA","SE","P"])
df.to_csv("${PRSICE_INPUT}", sep="\t", index=False)
print(f"  {len(df)} SNPs written")
PYEOF
else
  echo "  Sumstat already prepared — skipping"
fi

# ── Step 2: LD-clump at GWS (p < 5e-8) ──
CLUMP_OUT="${OUT_DIR}/${TRAIT}_clump_gws"
if [[ ! -f "${CLUMP_OUT}.clumped" ]]; then
  echo "  Clumping at p < 5e-8 (r2 < 0.001, 10Mb)..."
  "${PLINK1}" \
    --bfile "${UKB_GENO_STEP2}" \
    --clump "${PRSICE_INPUT}" \
    --clump-snp-field SNP \
    --clump-field    P \
    --clump-p1       5e-8 \
    --clump-p2       1 \
    --clump-r2       0.001 \
    --clump-kb       10000 \
    --out            "${CLUMP_OUT}" \
    --threads        "${N_THREADS}" 2>&1 | tail -5 || true
else
  echo "  Clump file already exists — skipping"
fi

if [[ ! -f "${CLUMP_OUT}.clumped" ]]; then
  echo "  WARNING: No clumped file produced — skipping scoring"
  exit 0
fi

N_CLUMP=$(tail -n +2 "${CLUMP_OUT}.clumped" | wc -l)
echo "  ${N_CLUMP} clumped SNPs (p < 5e-8)"

if [[ "${N_CLUMP}" -lt 1 ]]; then
  echo "  No genome-wide-significant SNPs — exiting"
  exit 0
fi

# ── Step 3: Build score file (SNP A1 BETA) ──
SCORE_FILE="${OUT_DIR}/${TRAIT}_score.txt"
if [[ ! -f "${SCORE_FILE}" ]]; then
  python3 - <<PYEOF
import pandas as pd
clump = pd.read_csv("${CLUMP_OUT}.clumped", sep=r"\s+")
sumstat = pd.read_csv("${PRSICE_INPUT}", sep="\t")
merged = clump[["SNP"]].merge(sumstat[["SNP","A1","BETA"]], on="SNP", how="left").dropna()
merged.to_csv("${SCORE_FILE}", sep="\t", index=False, header=False)
print(f"  Score file: {len(merged)} SNPs")
PYEOF
fi

# ── Step 4: Compute PRS per individual ──
PRS_OUT="${OUT_DIR}/${TRAIT}_prs"
if [[ ! -f "${PRS_OUT}.profile" ]]; then
  echo "  Computing PRS for all UKB individuals..."
  "${PLINK1}" \
    --bfile  "${UKB_GENO_STEP2}" \
    --score  "${SCORE_FILE}" 1 2 3 sum \
    --out    "${PRS_OUT}" \
    --threads "${N_THREADS}" 2>&1 | tail -5
else
  echo "  PRS profile already exists — skipping"
fi

echo ""
echo "[$(date)] ${TRAIT} PRS complete → ${OUT_DIR}/"
echo "  Score weights: ${SCORE_FILE}"
echo "  Profile:       ${PRS_OUT}.profile"
