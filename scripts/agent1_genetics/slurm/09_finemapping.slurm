#!/bin/bash
#------------------------------------------------------------
# SLURM: Fine-mapping with SuSiE at colocalized loci
# Runs after coloc results are combined
#------------------------------------------------------------
#SBATCH -J cad_finemapping
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 12:00:00
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

set +u  # conda activate scripts reference unbound vars (binutils ADDR2LINE)
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_R}"
set -u

mkdir -p "${AGENT1_OUT}/finemapping"

# Combine all coloc results into one file first
COLOC_DIR="${AGENT1_OUT}/coloc"
COMBINED_COLOC="${AGENT1_OUT}/coloc_results.csv"

echo "[$(date)] Combining coloc results..."
python - <<PYEOF
import os, pandas as pd

coloc_dir = "${COLOC_DIR}"
files = [f for f in os.listdir(coloc_dir) if f.endswith("_coloc.csv")]
dfs = [pd.read_csv(os.path.join(coloc_dir, f)) for f in files
       if os.path.getsize(os.path.join(coloc_dir, f)) > 0]
if dfs:
    combined = pd.concat(dfs, ignore_index=True)
    combined.to_csv("${COMBINED_COLOC}", index=False)
    n_coloc = combined['coloc_sig'].sum() if 'coloc_sig' in combined.columns else "unknown"
    print(f"Combined {len(dfs)} coloc files: {len(combined)} loci, {n_coloc} colocalized")
else:
    print("No coloc files found — check ${COLOC_DIR}")
PYEOF

echo "[$(date)] Running SuSiE fine-mapping..."
Rscript "${SCRIPTS_DIR}/agent1_genetics/R/run_finemapping.R" \
  --coloc_results "${COMBINED_COLOC}" \
  --sumstat_dir   "${AGENT1_OUT}/gwas_summary_stats/merged" \
  --ld_plink_dir  "${UKB_GENO_STEP2}" \
  --window_kb     500 \
  --pp4_thresh    0.8 \
  --out           "${AGENT1_OUT}/finemapping_credible_sets.csv"

# Copy key outputs to final results dir
cp "${COMBINED_COLOC}" "${AGENT1_OUT}/coloc_results.csv"
echo "[$(date)] Fine-mapping complete"
echo "Outputs: ${AGENT1_OUT}/finemapping_credible_sets.csv"
