#!/bin/bash
#------------------------------------------------------------
# Merge individual PRS profiles into one matrix
# Run after 10_prs.slurm (array 0-6) completes
#
# Submit: sbatch --dependency=afterok:<prs_jobid> 10b_prs_merge.slurm
#------------------------------------------------------------
#SBATCH -J cad_prs_merge
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH -A OTH26002
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=devanshpandey@utexas.edu

set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/work/07880/devansh/lonestar/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_PY}"
set -u

echo "[$(date)] Combining PRS profiles into summary matrix..."
python3 - <<PYEOF
import os, pandas as pd
from functools import reduce

prs_dir = "${AGENT1_OUT}/prs_weights"
traits = ["LDL_C", "HDL_C", "TRIGLY", "TOT_CHOL", "APOA1", "APOB", "LPA"]
dfs = []
for t in traits:
    f = os.path.join(prs_dir, t, f"{t}_prs.profile")
    if not os.path.exists(f):
        print(f"  Missing: {f}")
        continue
    df = pd.read_csv(f, sep=r"\s+")[["IID", "SCORESUM"]].rename(columns={"SCORESUM": t})
    dfs.append(df)
    print(f"  Loaded {t}: {len(df)} individuals")

if dfs:
    merged = reduce(lambda a, b: pd.merge(a, b, on="IID", how="inner"), dfs)
    out = os.path.join(prs_dir, "prs_all_traits.csv")
    merged.to_csv(out, index=False)
    print(f"\nCombined PRS matrix: {merged.shape[0]} individuals x {merged.shape[1]-1} traits")
    print(f"Written to: {out}")
else:
    print("ERROR: No PRS profiles found.")
    exit(1)
PYEOF

echo "[$(date)] PRS merge complete."
