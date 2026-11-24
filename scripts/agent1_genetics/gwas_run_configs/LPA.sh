# run_gwas.sh 1.2 — resolved configuration for LPA
# written 2026-09-04 14:15:00
# rerun a stage with:
#   /path/to/cad-genetics/scripts/run_gwas.sh --stage <step1|step2|merge> --config /path/to/cad-genetics/results/gwas/LPA/run_config.sh
SCRIPT_PATH="/path/to/cad-genetics/scripts/run_gwas.sh"
PHENO="/path/to/cad-genetics/results/agent1_genetics/phenotypes/pheno_quantitative.txt"
TRAIT="LPA"
OUTDIR="/path/to/cad-genetics/results/gwas/LPA"
COVAR="/path/to/cad-genetics/results/agent1_genetics/phenotypes/covariates.txt"
COVAR_COLS=""
CAT_COVAR=""
KEEP="/path/to/ukb/sample_eids.txt"
GENO_STEP1="/path/to/ukb_genotypes/merged_maf0.01_biallel_bbf_400k/merged_sub_chrom_maf0.01"
GENO_STEP2="/path/to/ukb_genotypes/merged_maf0.001_biallel_bbf_400k/merged_sub_chrom_maf0.001"
SAMPLE_FILE=""
PRUNED_SNPS="/path/to/ukb/GWAS_INPUT_FILES/prune_list_maf0.01.prune.in"
CHROMS="1 22"
MAF=""
MIN_MAC="20"
APPLY_RINT=1
GW_SIG="5e-8"
CLUMP_KB=1000
BSIZE_STEP1=1000
BSIZE_STEP2=400
ACCOUNT="YOUR_ALLOCATION"
PARTITION="normal"
TIME_STEP2="06:00:00"
NTASK_STEP2=16
THREADS=""
CONDA_ENV="regenie_env"
REGENIE_BIN=""
PLINK2_BIN=""
SCRATCH_ROOT="/path/to/scratch"
OUT_ROOT="/path/to/cad-genetics"
