#!/bin/bash
#SBATCH -J mr_robust
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 1:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
set -uo pipefail
_REPO="/path/to/cad-genetics"
source "${_REPO}/config/tacc_paths.sh"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_R}"
Rscript "${_REPO}/results/agent1_genetics/strengthen/mr_robustness.R"
echo "=== RESULT ==="; cat "${_REPO}/results/agent1_genetics/strengthen/mr_robustness/mr_robustness_CAD.csv"
