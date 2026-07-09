#!/bin/bash
#------------------------------------------------------------
# SLURM: Merge per-chromosome GWAS results and munge for LDSC
# Runs after jobs 03 and 04 complete.
#------------------------------------------------------------
#SBATCH -J cad_merge_munge
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 04:00:00
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

# This script uses two envs:
#   pyt_env  — for Python merging of REGENIE per-chr files
#   ldsc     — for munge_sumstats.py (activated per-section below)
set +u  # conda activate scripts reference unbound vars (binutils ADDR2LINE)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_PY}"
set -u

mkdir -p "${AGENT1_OUT}/gwas_summary_stats/merged" "${AGENT1_OUT}/ldsc_munged"

echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME:-shell}/${SLURM_JOB_ID:-0}] Merging per-chromosome GWAS results..."

merge_trait() {
  # $1 = full label, e.g. "lipids_LDL_C" or "outcomes_CAD"
  # $2 = directory containing per-chromosome files
  # REGENIE names files: {prefix}_chr{N}_{TRAIT}.regenie.gz
  #   e.g. lipids_chr1_LDL_C.regenie.gz
  local full_trait="$1"
  local trait_dir="$2"
  local out_prefix="${AGENT1_OUT}/gwas_summary_stats/merged/${full_trait}"

  # Split "lipids_LDL_C" → prefix="lipids", trait_name="LDL_C"
  local prefix="${full_trait%%_*}"       # everything before first _
  local trait_name="${full_trait#*_}"    # everything after first _

  # Skip if already merged
  if [[ -f "${out_prefix}.txt.gz" ]]; then
    echo "  Skipping ${full_trait} — merged file already exists"
    return 0
  fi

  echo "  Merging ${full_trait} (files: ${prefix}_chr*_${trait_name}.regenie.gz)..."

  # Verify chr1 exists before starting
  local chr1_file="${trait_dir}/${prefix}_chr1_${trait_name}.regenie.gz"
  if [[ ! -f "${chr1_file}" ]]; then
    echo "  ERROR: chr1 file not found: ${chr1_file}"
    echo "  Files present in ${trait_dir}:"; ls "${trait_dir}" | head -5
    return 1
  fi

  # REGENIE output columns: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P
  # Extract header from chr1 — use set +o pipefail because head -1 exits after one line,
  # causing gunzip to receive SIGPIPE (exit 141), which pipefail would treat as an error.
  set +o pipefail
  gunzip -c "${chr1_file}" | head -1 > "${out_prefix}.txt"
  set -o pipefail
  for chr in $(seq 1 22); do
    f="${trait_dir}/${prefix}_chr${chr}_${trait_name}.regenie.gz"
    if [[ -f "${f}" ]]; then
      gunzip -c "${f}" | tail -n +2 >> "${out_prefix}.txt"
    else
      echo "  WARNING: missing ${f}"
    fi
  done

  # Compress — prefer bgzip (enables tabix indexing); fall back to gzip
  # Use command -v so set -e doesn't fire on "command not found"
  if command -v bgzip &>/dev/null; then
    bgzip -@ 8 "${out_prefix}.txt"
    tabix -S 1 -s 1 -b 2 -e 2 "${out_prefix}.txt.gz" 2>/dev/null || true
  else
    gzip -f "${out_prefix}.txt"
    echo "  NOTE: bgzip not found — used gzip (tabix index skipped)"
  fi
  echo "  ${full_trait} merged: ${out_prefix}.txt.gz"
}

# Merge lipid traits
LIPID_TRAITS=(LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB LPA)
for trait in "${LIPID_TRAITS[@]}"; do
  merge_trait "lipids_${trait}" "${AGENT1_OUT}/gwas_summary_stats/lipids"
done

# Merge outcome traits — v2: added AF, PAD
OUTCOME_TRAITS=(CAD MI REVASC STROKE HF CV_DEATH MACE AF PAD)
for trait in "${OUTCOME_TRAITS[@]}"; do
  merge_trait "outcomes_${trait}" "${AGENT1_OUT}/gwas_summary_stats/outcomes"
done

echo "[$(date)] Munging sumstats for LDSC..."

# Step 1: Convert REGENIE LOG10P → Z-score using pyt_env (already active)
ALL_TRAITS=("${LIPID_TRAITS[@]/#/lipids_}" "${OUTCOME_TRAITS[@]/#/outcomes_}")
# Note: AF and PAD are included in OUTCOME_TRAITS above — no additional change needed
for trait in "${ALL_TRAITS[@]}"; do
  f="${AGENT1_OUT}/gwas_summary_stats/merged/${trait}.txt.gz"
  if [[ ! -f "${f}" ]]; then continue; fi

  [[ -f "${AGENT1_OUT}/ldsc_munged/${trait}_for_ldsc.txt" ]] && { echo "  Skipping ${trait} Z-score — already done"; continue; }
  echo "  Converting ${trait} to LDSC format..."
  python - <<PYEOF
import pandas as pd

# REGENIE step-2 output is space-separated (not tab).
df = pd.read_csv("${f}", sep=r"\s+", na_values=["NA"])
# munge_sumstats v1.0.1 requires a P column even with --signed-sumstats.
# Recover P exactly from REGENIE's LOG10P output.
df["Z"] = df["BETA"] / df["SE"]
df["P"] = 10.0 ** (-df["LOG10P"])
df = df.rename(columns={"ID": "SNP", "ALLELE1": "A1", "ALLELE0": "A2"})

out = "${AGENT1_OUT}/ldsc_munged/${trait}_for_ldsc.txt"
df[["SNP", "A1", "A2", "Z", "P", "N"]].to_csv(out, sep="\t", index=False)
print(f"  Written: {out}")
PYEOF
done

# Step 2: Switch to ldsc env for munge_sumstats.py
echo "[$(date)] Switching to ldsc env for munge_sumstats.py..."
set +u
conda activate "${CONDA_ENV_LDSC}"
set -u

for trait in "${ALL_TRAITS[@]}"; do
  input="${AGENT1_OUT}/ldsc_munged/${trait}_for_ldsc.txt"
  if [[ ! -f "${input}" ]]; then continue; fi
  [[ -f "${AGENT1_OUT}/ldsc_munged/${trait}.sumstats.gz" ]] && { echo "  Skipping ${trait} munge — already done"; continue; }
  echo "  Munging ${trait}..."
  MERGE_FLAG=""
  [[ -f "${LDSC_HAPMAP3_SNPS}" ]] && MERGE_FLAG="--merge-alleles ${LDSC_HAPMAP3_SNPS}"
  python "${MUNGE_PY}" \
    --sumstats "${input}" \
    ${MERGE_FLAG} \
    --out "${AGENT1_OUT}/ldsc_munged/${trait}" \
    --maf 0.01 \
    --chunksize 500000
  echo "  ${trait} munged"
done

echo "[$(date)] Merge and munge complete"
