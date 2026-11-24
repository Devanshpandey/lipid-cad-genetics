#!/bin/bash
#------------------------------------------------------------
# Statin-naive LDL-C / ApoB / triglyceride GWAS — REGENIE Step 1.
#
# WHY. The manuscript corrects LDL-C for lipid-lowering therapy (divide by 0.7)
# but does not correct ApoB, because no equivalent correction factor is
# established. A reviewer asked whether that asymmetry, rather than genuine
# collinearity, is what makes LDL-C and ApoB non-identifiable in the joint model.
# The clean test is to drop the correction entirely and analyse only participants
# who report no lipid-lowering medication at baseline: there, both traits are
# measured on the same untreated footing. If the conditional F-statistics still
# collapse, non-identifiability is a property of the genetics, not of the
# medication handling.
#------------------------------------------------------------
#SBATCH -J snaive1
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 64
#SBATCH -t 24:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/snaive1_%j.out
#SBATCH -e /path/to/cad-genetics/snaive1_%j.err
set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

PH="${AGENT1_OUT}/phenotypes"
KEEP="${PH}/regenie_keep_statinnaive.txt"
RAW=/path/to/ukb_phenotypes

echo "[$(date)] building the statin-naive keep list"
conda activate "${CONDA_ENV_R}"
PH_DIR="${PH}" Rscript - <<'RS'
library(data.table)
raw <- "/path/to/ukb_phenotypes"
ph  <- Sys.getenv("PH_DIR")
rd <- function(fid) {
  d <- fread(file.path(raw, sprintf("fid%s.csv", fid)), na.strings = c("", "NA"))
  col <- grep(sprintf("^\"?%s-0\\.0\"?$", fid), names(d), value = TRUE)[1]
  setnames(d, c("eid", col), c("IID", "v")); d[, .(IID = as.character(IID), v)]
}
# Medication fields are sex-split: 6153 asked of women, 6177 of men. Reading only
# one of them is the ascertainment bug this project already fixed once.
st <- merge(rd("6153")[, .(IID, f = v)], rd("6177")[, .(IID, m = v)], by = "IID", all = TRUE)
st[, statin := as.integer((!is.na(f) & f == 1) | (!is.na(m) & m == 1))]
st[is.na(f) & is.na(m), statin := NA_integer_]
keep <- fread(file.path(ph, "regenie_keep.txt"), header = FALSE, col.names = c("FID","IID"))
keep[, IID := as.character(IID)]
naive <- merge(keep, st[, .(IID, statin)], by = "IID")[statin == 0]
cat(sprintf("[naive] cohort %d -> statin-naive %d (%.1f%%); statin status missing for %d\n",
            nrow(keep), nrow(naive), 100*nrow(naive)/nrow(keep),
            nrow(keep) - nrow(merge(keep, st[!is.na(statin), .(IID)], by = "IID"))))
fwrite(naive[, .(FID, IID)], file.path(ph, "regenie_keep_statinnaive.txt"),
       sep = " ", col.names = FALSE)
RS

echo "[$(date)] keep list: $(wc -l < ${KEEP}) participants"

conda activate "${CONDA_ENV_REGENIE}"
mkdir -p "${TACC_SCRATCH}/regenie_tmp" "${AGENT1_OUT}/regenie_step1_statinnaive"
regenie \
  --step 1 \
  --bed "${UKB_GENO_STEP1}" \
  --extract "${PH}/qc_pruned_snps.txt" \
  --keep "${KEEP}" \
  --phenoFile "${PH}/pheno_quantitative.txt" \
  --phenoColList LDL_C,APOB,TRIGLY \
  --covarFile "${PH}/covariates.txt" \
  --bsize 1000 \
  --lowmem \
  --lowmem-prefix "${TACC_SCRATCH}/regenie_tmp/step1_naive" \
  --threads "${SLURM_NTASKS}" \
  --gz \
  --out "${AGENT1_OUT}/regenie_step1_statinnaive/step1_naive"
echo "[$(date)] step 1 complete"
