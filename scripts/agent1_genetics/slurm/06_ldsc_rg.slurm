#!/bin/bash
#------------------------------------------------------------
# SLURM: LDSC genetic correlations
#   Block 1 — lipids × outcomes  (7 × 7 = 49 pairs)
#   Block 2 — lipids × lipids    (7 × 7 = 49 pairs, upper triangle = 28 unique)
# Total: 98 LDSC runs → two output CSVs + one combined matrix
#------------------------------------------------------------
#SBATCH -J cad_ldsc_rg
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 06:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail
# Resolve repo root: SLURM_SUBMIT_DIR when run via sbatch, hardcoded path as fallback
# (BASH_SOURCE[0] points to /var/spool in SLURM jobs — cannot use dirname approach)
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

set +u  # conda activate scripts reference unbound vars (binutils ADDR2LINE)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_LDSC}"
set -u

OUT_DIR="${AGENT1_OUT}/ldsc_rg"
mkdir -p "${OUT_DIR}"

echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME:-shell}/${SLURM_JOB_ID:-0}] Running LDSC genetic correlations..."

MUNGE_DIR="${AGENT1_OUT}/ldsc_munged"
LIPID_TRAITS=(lipids_LDL_C lipids_HDL_C lipids_TRIGLY lipids_TOT_CHOL lipids_APOA1 lipids_APOB lipids_LPA)
# v2: added outcomes_AF (atrial fibrillation) and outcomes_PAD (peripheral artery disease)
OUTCOME_TRAITS=(outcomes_CAD outcomes_MI outcomes_REVASC outcomes_STROKE outcomes_HF outcomes_CV_DEATH outcomes_MACE outcomes_AF outcomes_PAD)

# Detect weights prefix once
WEIGHTS_PREFIX="${LDSC_WEIGHTS_DIR}/weights.hm3_noMHC."
if [[ ! -f "${LDSC_WEIGHTS_DIR}/weights.hm3_noMHC.1.l2.ldscore.gz" ]]; then
  WEIGHTS_PREFIX="${LDSC_REF_PREFIX}"
  echo "  NOTE: using baselineLD files for weights (standard weights not found)"
fi

# ── helper: run one LDSC rg call ────────────────────────────
run_rg() {
  local exposure="$1"
  local comparator_list="$2"   # comma-separated .sumstats.gz paths
  local out_prefix="$3"
  python "${LDSC_PY}" \
    --rg "${exposure},${comparator_list}" \
    --ref-ld-chr "${LDSC_REF_PREFIX}" \
    --w-ld-chr "${WEIGHTS_PREFIX}" \
    --out "${out_prefix}"
}

# ── Block 1: lipids × outcomes (7 × 7 = 49 pairs) ──────────
echo "[$(date)] === Block 1: lipids × outcomes ==="
for lipid in "${LIPID_TRAITS[@]}"; do
  exposure="${MUNGE_DIR}/${lipid}.sumstats.gz"
  [[ ! -f "${exposure}" ]] && { echo "WARNING: missing ${exposure}"; continue; }

  outcome_list=""
  for outcome in "${OUTCOME_TRAITS[@]}"; do
    f="${MUNGE_DIR}/${outcome}.sumstats.gz"
    [[ -f "${f}" ]] && outcome_list="${outcome_list},${f}"
  done
  outcome_list="${outcome_list#,}"
  [[ -z "${outcome_list}" ]] && { echo "WARNING: no outcome sumstats found"; continue; }

  echo "[$(date)]   ${lipid} vs all outcomes..."
  run_rg "${exposure}" "${outcome_list}" "${OUT_DIR}/${lipid}_vs_outcomes"
done

# ── Block 2: lipids × lipids (upper triangle, 21 unique pairs) ──
echo "[$(date)] === Block 2: lipids × lipids ==="
for i in "${!LIPID_TRAITS[@]}"; do
  lipid="${LIPID_TRAITS[$i]}"
  exposure="${MUNGE_DIR}/${lipid}.sumstats.gz"
  [[ ! -f "${exposure}" ]] && { echo "WARNING: missing ${exposure}"; continue; }

  # Build list of all OTHER lipid traits (full 7, LDSC deduplicates self-correlation)
  other_list=""
  for other in "${LIPID_TRAITS[@]}"; do
    f="${MUNGE_DIR}/${other}.sumstats.gz"
    [[ -f "${f}" && "${other}" != "${lipid}" ]] && other_list="${other_list},${f}"
  done
  other_list="${other_list#,}"
  [[ -z "${other_list}" ]] && continue

  echo "[$(date)]   ${lipid} vs all other lipids..."
  run_rg "${exposure}" "${other_list}" "${OUT_DIR}/${lipid}_vs_lipids"
done

# ── Compile all results into CSV matrices ───────────────────
# Switch back to pyt_env (Python 3 + polars) — ldsc env is Python 2
echo "[$(date)] Compiling rg matrices..."
set +u
conda activate "${CONDA_ENV_PY}"
set -u

python - <<'PYEOF'
# -*- coding: utf-8 -*-
import os, re, polars as pl

out_dir  = os.environ["AGENT1_OUT"] + "/ldsc_rg"
agent1   = os.environ["AGENT1_OUT"]

def parse_ldsc_log(log_path, exposure_label):
    rows = []
    try:
        with open(log_path) as fh:
            content = fh.read()
    except FileNotFoundError:
        return rows
    # LDSC prints a summary table; grab each result line after the header
    for block in re.finditer(
        r"p1\s+p2\s+rg\s+se\s+z\s+p[^\n]*\n((?:[^\n]+\n?)*)",
        content
    ):
        for line in block.group(1).strip().split("\n"):
            parts = line.split()
            if len(parts) >= 6 and not parts[0].startswith("p1"):
                rows.append({
                    "exposure": exposure_label,
                    "outcome":  os.path.basename(parts[1]).replace(".sumstats.gz",""),
                    "rg":  float(parts[2]) if parts[2] not in ("nan","NA") else None,
                    "se":  float(parts[3]) if parts[3] not in ("nan","NA") else None,
                    "z":   float(parts[4]) if parts[4] not in ("nan","NA") else None,
                    "p":   float(parts[5]) if parts[5] not in ("nan","NA") else None,
                })
    return rows

lipid_traits   = ["lipids_LDL_C","lipids_HDL_C","lipids_TRIGLY","lipids_TOT_CHOL",
                  "lipids_APOA1","lipids_APOB","lipids_LPA"]
outcome_traits = ["outcomes_CAD","outcomes_MI","outcomes_REVASC","outcomes_STROKE",
                  "outcomes_HF","outcomes_CV_DEATH","outcomes_MACE",
                  "outcomes_AF","outcomes_PAD"]

# lipids × outcomes matrix
rows_lo = []
for lipid in lipid_traits:
    log = os.path.join(out_dir, f"{lipid}_vs_outcomes.log")
    rows_lo.extend(parse_ldsc_log(log, lipid))

# lipids × lipids matrix
rows_ll = []
for lipid in lipid_traits:
    log = os.path.join(out_dir, f"{lipid}_vs_lipids.log")
    rows_ll.extend(parse_ldsc_log(log, lipid))

# Write separate CSVs
for rows, fname, label in [
    (rows_lo, "ldsc_rg_lipids_vs_outcomes.csv", "lipids×outcomes"),
    (rows_ll, "ldsc_rg_lipids_vs_lipids.csv",   "lipids×lipids"),
]:
    if rows:
        df = pl.DataFrame(rows)
        path = os.path.join(agent1, fname)
        df.write_csv(path)
        print(f"  {label}: {path} ({len(df)} pairs)")
    else:
        print(f"  WARNING: no results parsed for {label}")

# Combined matrix (all pairs)
all_rows = rows_lo + rows_ll
if all_rows:
    df_all = pl.DataFrame(all_rows)
    combined = os.path.join(agent1, "ldsc_rg_matrix.csv")
    df_all.write_csv(combined)
    print(f"  Combined: {combined} ({len(df_all)} total pairs)")
PYEOF

echo "[$(date)] LDSC rg complete."
echo "  Outputs:"
echo "    ${AGENT1_OUT}/ldsc_rg_lipids_vs_outcomes.csv"
echo "    ${AGENT1_OUT}/ldsc_rg_lipids_vs_lipids.csv"
echo "    ${AGENT1_OUT}/ldsc_rg_matrix.csv  (combined)"
