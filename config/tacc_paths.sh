#!/usr/bin/env bash
# ============================================================
# TACC Lonestar6 — path configuration
# Source this file at the top of every SLURM script.
# No edits needed — all paths confirmed from your environment.
# ============================================================

# Guard against conda.sh "HOST: unbound variable" under set -u
# (conda 4.x activation script references $HOST; ensure it's set before sourcing)
export HOST="${HOST:-${HOSTNAME:-$(hostname)}}"

# --- Your TACC identity ---
TACC_USER="devansh"
TACC_HOST="ls6.tacc.utexas.edu"
TACC_ALLOCATION="OTH26002"

# --- Project directories ---
TACC_WORK="/work/07880/devansh/lonestar/cad-genetics"
TACC_SCRATCH="/scratch/07880/devansh/cad-genetics"

# ---------------------------------------------------------------
# EXISTING CONDA ENVIRONMENTS  (do not recreate these)
# ---------------------------------------------------------------
CONDA_BASE="/work/07880/devansh/anaconda3"
CONDA_ENV_REGENIE="regenie_env"   # REGENIE + dependencies
CONDA_ENV_R="r_env"               # R + stats packages
CONDA_ENV_LDSC="ldsc"             # LDSC Python environment
CONDA_ENV_PY="pyt_env"            # General Python (pandas, polars, etc.)

# ---------------------------------------------------------------
# EXISTING TOOL BINARIES  (absolute paths — no module load needed)
# ---------------------------------------------------------------
PLINK2="/work/07880/devansh/lonestar/software/plink2"
BCFTOOLS="/work/07880/devansh/lonestar/software/bcftools/bcftools"

# LDSC scripts
LDSC_PY="/work/07880/devansh/lonestar/software/ldsc/ldsc.py"
MUNGE_PY="/work/07880/devansh/lonestar/software/ldsc/munge_sumstats.py"

# ---------------------------------------------------------------
# UKB GENOTYPE DATA  (corral — read-only)
# ---------------------------------------------------------------
UKB_GENO_ROOT="/corral/utexas/UKB-Imaging-Genetics/UKB_GENOTYPE_QC_400k"

# Single merged BED/BIM/FAM — MAF > 0.001 — used for REGENIE Step 2 (discovery)
UKB_GENO_STEP2="${UKB_GENO_ROOT}/merged_maf0.001_biallel_bbf_400k/merged_sub_chrom_maf0.001"

# MAF > 0.01 merged file — used for REGENIE Step 1 null model (LD pruning)
UKB_GENO_STEP1="${UKB_GENO_ROOT}/merged_maf0.01_biallel_bbf_400k/merged_sub_chrom_maf0.01"

# QC-passed white British EID list (one EID per line, no header)
UKB_SAMPLE_EID="${UKB_GENO_ROOT}/geno_qc_eids_400k_white_british_2023-11-29.txt"

# ---------------------------------------------------------------
# UKB WES DATA  (corral-repl — read-only, all-chr merged pgen)
# Confirmed path: /corral-repl/utexas/UKB-Imaging-Genetics/UKB_GENOTYPE/ukb_pfile/
#   wes_allchr.pgen  (40 GB)   — WES genotypes, all chromosomes
#   wes_allchr.pvar  (863 MB)  — variant info
#   wes_allchr.psam  (8.5 MB)  — sample info
#   joyce_annotated_vcf/       — VEP-annotated variant files (used for burden masks)
#   helper_files/              — supplementary files
# ---------------------------------------------------------------
UKB_WES_ROOT="/corral-repl/utexas/UKB-Imaging-Genetics/UKB_GENOTYPE/ukb_pfile"
UKB_WES_PFILE="${UKB_WES_ROOT}/wes_allchr"          # prefix for --pfile in REGENIE/PLINK2
UKB_WES_JOYCE_VCF="${UKB_WES_ROOT}/joyce_annotated_vcf"  # VEP annotation VCFs

# ---------------------------------------------------------------
# UKB PHENOTYPE DATA  (corral — read-only)
# Each field → fid{field_id}.csv  |  columns: eid, {field}-0.0, {field}-1.0
# ---------------------------------------------------------------
UKB_PHENO_DIR="/corral/utexas/UKB-Imaging-Genetics/temp_imaging_data/pheno_split_into_files_011924"
UKB_ICD_FILE="${UKB_PHENO_DIR}/binary_ICD_011924.txt"

# ---------------------------------------------------------------
# LDSC REFERENCE  (pre-existing at software install path)
# baselineLD v2.2 files: baselineLD.{chr}.l2.ldscore.gz
# Use for both --ref-ld-chr and --w-ld-chr in ldsc.py --rg
# ---------------------------------------------------------------
LDSC_REF_DIR="/work/07880/devansh/lonestar/software/ldsc/data/ldsc_references"
# Full prefix passed to --ref-ld-chr / --w-ld-chr (chr number appended by ldsc.py)
LDSC_REF_PREFIX="${LDSC_REF_DIR}/baselineLD."
# Weights dir (1000G HapMap3 no-MHC) — ls its contents to confirm filename prefix
LDSC_WEIGHTS_DIR="${LDSC_REF_DIR}/1000G_Phase3_weights_hm3_no_MHC"
# HapMap3 SNP list for munge_sumstats --merge-alleles (optional; skipped if absent)
LDSC_HAPMAP3_SNPS="${LDSC_REF_DIR}/w_hm3.snplist"

# ---------------------------------------------------------------
# DERIVED PATHS  (do not edit)
# ---------------------------------------------------------------
export SCRIPTS_DIR="${TACC_WORK}/scripts"
export RESULTS_DIR="${TACC_WORK}/results"
export AGENT1_OUT="${RESULTS_DIR}/agent1_genetics"
export AGENT2_OUT="${RESULTS_DIR}/agent2_genes"
export AGENT3_OUT="${RESULTS_DIR}/agent3_networks"
export AGENT4_OUT="${RESULTS_DIR}/agent4_subtypes"
export LOGS_DIR="${TACC_WORK}/logs"
export DATA_DIR="${TACC_WORK}/data"

# Burden mask annotation files — defined here so AGENT1_OUT is already set
# (produced by prepare_burden_masks.slurm + merge_burden_masks.slurm)
export UKB_BURDEN_ANNO_DIR="${AGENT1_OUT}/burden_masks"
export UKB_BURDEN_ANNO_FILE="${UKB_BURDEN_ANNO_DIR}/ukb_wes_anno.txt"
export UKB_BURDEN_SETLIST="${UKB_BURDEN_ANNO_DIR}/ukb_wes_setlist.txt"
export UKB_BURDEN_MASKS="${UKB_BURDEN_ANNO_DIR}/ukb_wes_masks.txt"

# ---------------------------------------------------------------
# HELPER: activate a conda env in a SLURM script
# Usage:  conda_activate "${CONDA_ENV_R}"
# ---------------------------------------------------------------
conda_activate() {
  source "${CONDA_BASE}/etc/profile.d/conda.sh"
  conda activate "$1"
}
