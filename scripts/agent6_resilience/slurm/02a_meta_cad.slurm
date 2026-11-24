#!/bin/bash
#------------------------------------------------------------
# Agent 6 - fixed-effect meta-analysis of UKB-free CAD GWAS (MVP EUR + Nikpay 2015 + FinnGen R11)
# Output: results/agent6_resilience/prs/cad_ukbfree_meta_sst.txt (+ .summary.json)
# Submit FROM THE REPO ROOT:  sbatch scripts/agent6_resilience/slurm/02a_meta_cad.slurm
#------------------------------------------------------------
#SBATCH -J a6_meta_cad
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o logs/a6_meta_cad_%j.out
#SBATCH -e logs/a6_meta_cad_%j.err
set -euo pipefail
_REPO="/path/to/cad-genetics"
source "${_REPO}/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_PY}"; set -u
OUT_DIR="${RESULTS_DIR}/agent6_resilience/prs"; mkdir -p "${OUT_DIR}"
EXT="${_REPO}/data/external"
python3 "${_REPO}/scripts/agent6_resilience/prs/meta_cad_ukbfree.py" \
  --mvp "${EXT}/mvp_cad/mvp_cad_EUR.txt.gz" \
  --nikpay "${EXT}/cardiogram/cad.add.160614.website.txt" \
  --finngen "${EXT}/finngen/finngen_R11_I9_CHD.gz" \
  --n_finngen_cases 51098 --n_finngen_controls 402635 \
  --out "${OUT_DIR}/cad_ukbfree_meta_sst.txt"
md5sum "${EXT}/mvp_cad/mvp_cad_EUR.txt.gz" "${EXT}/cardiogram/cad.add.160614.website.txt" "${EXT}/finngen/finngen_R11_I9_CHD.gz" \
  "${OUT_DIR}/cad_ukbfree_meta_sst.txt" > "${OUT_DIR}/cad_ukbfree_meta_sst.md5"
cat "${OUT_DIR}/cad_ukbfree_meta_sst.txt.summary.json"
