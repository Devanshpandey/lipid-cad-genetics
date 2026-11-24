#!/bin/bash
#------------------------------------------------------------
# Agent 6 - one-time extraction of all PRS-CS-weighted SNPs (union of every
# *_all_chr.txt weight file) from the 1.5 TB corral bed into a small subset bed
# on scratch. The full-file read costs ~4 h once; every scoring job afterwards
# reads ~120 GB and finishes in minutes. Scoring jobs use this subset when present.
# Submit FROM THE REPO ROOT:  sbatch scripts/agent6_resilience/slurm/02c_extract_prs_bed.slurm
#------------------------------------------------------------
#SBATCH -J a6_extract_bed
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 08:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o logs/a6_extract_bed_%j.out
#SBATCH -e logs/a6_extract_bed_%j.err
set -euo pipefail
_REPO="/path/to/cad-genetics"
source "${_REPO}/config/tacc_paths.sh"
PRS_DIR="${RESULTS_DIR}/agent6_resilience/prs"
SUB_DIR="${TACC_SCRATCH}/agent6_prs_subset"; mkdir -p "${SUB_DIR}"
SNPS="${SUB_DIR}/prs_snps_union.txt"
ls "${PRS_DIR}"/*_all_chr.txt
awk '{print $2}' "${PRS_DIR}"/*_all_chr.txt | sort -u > "${SNPS}"
echo "union of weighted SNPs: $(wc -l < "${SNPS}")"
"${PLINK2}" --bfile "${UKB_GENO_STEP2}" --keep "${UKB_SAMPLE_EID}" --extract "${SNPS}" \
  --make-bed --threads 16 --memory 200000 --out "${SUB_DIR}/ukb400k_prs_snps"
wc -l "${SUB_DIR}/ukb400k_prs_snps.bim" "${SUB_DIR}/ukb400k_prs_snps.fam"
md5sum "${SNPS}" "${SUB_DIR}/ukb400k_prs_snps.bim" > "${SUB_DIR}/ukb400k_prs_snps.md5"
echo "[$(date)] subset bed ready: ${SUB_DIR}/ukb400k_prs_snps"
