#!/bin/bash
#------------------------------------------------------------
# Agent 6 - CAD PRS from the UKB-free meta-analysis (MVP EUR + Nikpay 2015 + FinnGen R11), PRS-CS-auto.
# Array over 22 chromosomes; needs 02a_meta_cad.slurm first. Follow with 02b_prscs_cad_score.slurm.
# Set SST_TAG=nikpay2015 to build the Nikpay-only sensitivity PRS instead.
# Submit FROM THE REPO ROOT:  sbatch --array=1-22 scripts/agent6_resilience/slurm/02_prscs_cad.slurm
#------------------------------------------------------------
#SBATCH -J a6_prscs_cad
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 04:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o logs/a6_prscs_cad_%A_%a.out
#SBATCH -e logs/a6_prscs_cad_%A_%a.err
set -euo pipefail
_REPO="/path/to/cad-genetics"
source "${_REPO}/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_PY}"; set -u
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PRSCS_DIR="/path/to/software/PRScs"
CHR="${SLURM_ARRAY_TASK_ID}"
OUT_DIR="${RESULTS_DIR}/agent6_resilience/prs"; mkdir -p "${OUT_DIR}"
SST_TAG="${SST_TAG:-ukbfree}"
if [[ "${SST_TAG}" == "ukbfree" ]]; then
  SST="${OUT_DIR}/cad_ukbfree_meta_sst.txt"
  [[ -s "${SST}" ]] || { echo "ERROR: run 02a_meta_cad.slurm first (${SST} missing)"; exit 1; }
  # PRS-CS uses one N: the median per-SNP effective N of the meta-analysis
  N_GWAS=$(python3 -c "import json;print(json.load(open('${SST}.summary.json'))['n_gwas_for_prscs'])")
  PREFIX="cad_ukbfree"
else
  SST="${OUT_DIR}/nikpay2015_prscs_sst.txt"; N_GWAS=162970; PREFIX="cad_nikpay"
  if [[ ! -s "${SST}" ]]; then
    python3 "${_REPO}/scripts/agent6_resilience/prs/prepare_nikpay_sst.py" "${_REPO}/data/external/cardiogram/cad.add.160614.website.txt" "${SST}" | tee "${OUT_DIR}/sst_prep.log"
  fi
fi
echo "sumstats=${SST} n_gwas=${N_GWAS} prefix=${PREFIX}"
EXPECTED="${OUT_DIR}/${PREFIX}_pst_eff_a1_b0.5_phiauto_chr${CHR}.txt"
if [[ -s "${EXPECTED}" ]]; then echo "exists: ${EXPECTED}"; exit 0; fi
python3 "${PRSCS_DIR}/PRScs.py" --ref_dir="${PRSCS_DIR}/ldblk_ukbb_eur" --bim_prefix="${UKB_GENO_STEP2}" \
  --sst_file="${SST}" --n_gwas="${N_GWAS}" --chrom="${CHR}" --out_dir="${OUT_DIR}/${PREFIX}" --seed=42
[[ -s "${EXPECTED}" ]] || { echo "ERROR: no output for chr${CHR}"; exit 1; }
echo "[$(date)] chr${CHR} done: $(wc -l < "${EXPECTED}") SNPs"
