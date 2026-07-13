#!/bin/bash
#------------------------------------------------------------
# SLURM: REGENIE Step 1 — fit null model (whole-genome regression)
# Runs once; outputs predictions used by all Step 2 jobs.
# Uses --lowmem to handle UKB scale (~500K samples, ~800K SNPs).
#------------------------------------------------------------
#SBATCH -J cad_regenie_step1
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 24:00:00
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

source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_REGENIE}"

mkdir -p "${TACC_SCRATCH}/regenie_tmp" "${AGENT1_OUT}/regenie_step1"

echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME:-shell}/${SLURM_JOB_ID:-0}] REGENIE Step 1 — null model fitting"
echo "Threads: ${SLURM_NTASKS}"

# Step 1 uses a pruned/LD-independent SNP set for computational efficiency.
# Typically: MAF > 0.01, genotyping rate > 0.99, HWE p > 1e-15, r2 < 0.1 (window 1000, step 100)
# Pre-generate this SNP list with PLINK2 if not already done.
PRUNE_SNP_LIST="${AGENT1_OUT}/phenotypes/qc_pruned_snps.txt"

if [[ ! -f "${PRUNE_SNP_LIST}" ]]; then
  echo "[$(date)] Generating pruned SNP list for Step 1 (from MAF 0.01 merged file)..."
  # UKB_GENO_STEP1 points to the merged MAF>0.01 BED file — already MAF filtered.
  # Additional LD pruning gives the ~100-500K SNP set REGENIE expects for Step 1.
  "${PLINK2}" \
    --bed "${UKB_GENO_STEP1}.bed" \
    --bim "${UKB_GENO_STEP1}.bim" \
    --fam "${UKB_GENO_STEP1}.fam" \
    --keep "${AGENT1_OUT}/phenotypes/regenie_keep.txt" \
    --maf 0.01 \
    --geno 0.01 \
    --hwe 1e-15 \
    --indep-pairwise 1000 100 0.1 \
    --out "${AGENT1_OUT}/phenotypes/qc_pruned" \
    --threads "${SLURM_NTASKS}"
  mv "${AGENT1_OUT}/phenotypes/qc_pruned.prune.in" "${PRUNE_SNP_LIST}"
  echo "[$(date)] Pruned SNP list: $(wc -l < ${PRUNE_SNP_LIST}) SNPs"
fi

# REGENIE Step 1 — quantitative traits (lipids)
# Uses merged MAF>0.01 BED file with LD-pruned SNP set
echo "[$(date)] Running REGENIE Step 1 (quantitative traits)..."
regenie \
  --step 1 \
  --bed "${UKB_GENO_STEP1}" \
  --extract "${PRUNE_SNP_LIST}" \
  --keep "${AGENT1_OUT}/phenotypes/regenie_keep.txt" \
  --phenoFile "${AGENT1_OUT}/phenotypes/pheno_quantitative.txt" \
  --covarFile "${AGENT1_OUT}/phenotypes/covariates.txt" \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix "${TACC_SCRATCH}/regenie_tmp/step1_qt" \
  --threads "${SLURM_NTASKS}" \
  --gz \
  --out "${AGENT1_OUT}/regenie_step1/step1_qt"

# REGENIE Step 1 — binary traits (CAD/MACE outcomes)
echo "[$(date)] Running REGENIE Step 1 (binary traits)..."
regenie \
  --step 1 \
  --bed "${UKB_GENO_STEP1}" \
  --extract "${PRUNE_SNP_LIST}" \
  --keep "${AGENT1_OUT}/phenotypes/regenie_keep.txt" \
  --phenoFile "${AGENT1_OUT}/phenotypes/pheno_binary.txt" \
  --covarFile "${AGENT1_OUT}/phenotypes/covariates.txt" \
  --bsize 1000 \
  --bt \
  --lowmem \
  --lowmem-prefix "${TACC_SCRATCH}/regenie_tmp/step1_bt" \
  --threads "${SLURM_NTASKS}" \
  --gz \
  --out "${AGENT1_OUT}/regenie_step1/step1_bt"

echo "[$(date)] REGENIE Step 1 complete"
echo "Prediction files in: ${AGENT1_OUT}/regenie_step1/"
