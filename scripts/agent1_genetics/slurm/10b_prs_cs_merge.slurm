#!/bin/bash
#------------------------------------------------------------
# Merge PRS-CS per-chr posterior weights and score all UKB individuals
# Run after 10_prs_cs.slurm (array 0-153) completes.
#
# For each trait:
#   1. Cat per-chr posterior weights into one file
#   2. Score all UKB individuals with plink1.9 --score
#   3. Combine all 7 trait profiles into prs_all_traits.csv
#
# Submit: sbatch --dependency=afterok:<prscs_jobid> 10b_prs_cs_merge.slurm
#------------------------------------------------------------
#SBATCH -J cad_prscs_merge
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
_REPO="${SLURM_SUBMIT_DIR:-/work/07880/devansh/lonestar/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

PLINK1="/work/07880/devansh/lonestar/software/plink"
[[ -f "${PLINK1}" ]] || { echo "ERROR: plink1.9 not found"; exit 1; }

set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_PY}"
set -u

N_THREADS="${SLURM_NTASKS:-16}"
PRSCS_DIR="${AGENT1_OUT}/prs_weights_prscs"
PRS_OUT="${AGENT1_OUT}/prs_weights"
mkdir -p "${PRS_OUT}"

LIPID_TRAITS=(LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB LPA)

echo "[$(date)] === Merging PRS-CS weights and scoring ==="

for TRAIT in "${LIPID_TRAITS[@]}"; do
  echo ""
  echo "[$(date)] ---- ${TRAIT} ----"
  TRAIT_DIR="${PRSCS_DIR}/${TRAIT}"

  # Check all 22 chr files exist
  MISSING=0
  for CHR in $(seq 1 22); do
    f="${TRAIT_DIR}/${TRAIT}_prscs_pst_eff_a1_b0.5_phiauto_chr${CHR}.txt"
    [[ -f "$f" ]] || { echo "  WARNING: missing chr${CHR}"; MISSING=$((MISSING+1)); }
  done
  if [[ "${MISSING}" -gt 0 ]]; then
    echo "  SKIP: ${MISSING} chromosome(s) missing for ${TRAIT}"
    continue
  fi

  # Merge per-chr weights
  MERGED_WEIGHTS="${TRAIT_DIR}/${TRAIT}_prscs_weights_all.txt"
  if [[ ! -f "${MERGED_WEIGHTS}" ]]; then
    echo "  Merging per-chr weights..."
    cat "${TRAIT_DIR}/${TRAIT}_prscs_pst_eff_a1_b0.5_phiauto_chr"*.txt > "${MERGED_WEIGHTS}"
    echo "  $(wc -l < "${MERGED_WEIGHTS}") total SNPs"
  fi

  # PRS-CS output format: CHR SNP BP A1 A2 POSTERIOR_BETA (cols 1-6)
  # plink --score needs: SNP A1 BETA → cols 2, 4, 6
  SCORE_FILE="${TRAIT_DIR}/${TRAIT}_prscs_score.txt"
  awk '{print $2, $4, $6}' "${MERGED_WEIGHTS}" > "${SCORE_FILE}"

  # Score with plink
  PRS_PROFILE="${PRS_OUT}/${TRAIT}/${TRAIT}_prscs_prs"
  mkdir -p "${PRS_OUT}/${TRAIT}"
  if [[ ! -f "${PRS_PROFILE}.profile" ]]; then
    echo "  Scoring all UKB individuals..."
    "${PLINK1}" \
      --bfile "${UKB_GENO_STEP2}" \
      --score "${SCORE_FILE}" 1 2 3 sum \
      --out "${PRS_PROFILE}" \
      --threads "${N_THREADS}" 2>&1 | tail -5
  else
    echo "  Profile already exists — skipping"
  fi

  echo "  Done: ${PRS_PROFILE}.profile"
done

# Combine all profiles
echo ""
echo "[$(date)] Combining PRS profiles into summary matrix..."
python3 - <<PYEOF
import os, pandas as pd
from functools import reduce

prs_dir = "${PRS_OUT}"
traits = ["LDL_C", "HDL_C", "TRIGLY", "TOT_CHOL", "APOA1", "APOB", "LPA"]
dfs = []
for t in traits:
    f = os.path.join(prs_dir, t, f"{t}_prscs_prs.profile")
    if not os.path.exists(f):
        print(f"  Missing: {f}")
        continue
    df = pd.read_csv(f, sep=r"\s+")[["IID", "SCORESUM"]].rename(columns={"SCORESUM": f"{t}_prscs"})
    dfs.append(df)
    print(f"  Loaded {t}: {len(df)} individuals")

if dfs:
    merged = reduce(lambda a, b: pd.merge(a, b, on="IID", how="inner"), dfs)
    out = os.path.join(prs_dir, "prs_all_traits_prscs.csv")
    merged.to_csv(out, index=False)
    print(f"\nCombined PRS-CS matrix: {merged.shape[0]} individuals x {merged.shape[1]-1} traits")
    print(f"Written to: {out}")
else:
    print("ERROR: No PRS-CS profiles found.")
    exit(1)
PYEOF

echo "[$(date)] === PRS-CS merge complete ==="
