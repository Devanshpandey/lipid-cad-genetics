#!/bin/bash
#------------------------------------------------------------
# Lp(a) GWAS with the too-low / too-high participants kept (field 30796), alongside a reproduction control.
#
# Sensitivity analysis only: nothing under merged/ or merged_statinfix/ is touched.
# Stage A (default): phenotype, step 1, chr 6, and the control gate.
#   sbatch 08ac_lpa_censored_gwas.slurm
# Stage B: the other chromosomes, merge, and the sumstat dirs 08ad consumes.
#   sbatch --dependency=afterok:<A> --export=ALL,STAGE=B 08ac_lpa_censored_gwas.slurm
# Step 2 mirrors 03_regenie_step2_lipids.slurm except that the scores are precomputed, so no --apply-rint.
#------------------------------------------------------------
#SBATCH -J lpacens_gwas
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/lpacens_gwas_%j.out
#SBATCH -e /path/to/cad-genetics/logs/lpacens_gwas_%j.err
set -euo pipefail
umask 077
_REPO=/path/to/cad-genetics
source "${_REPO}/config/tacc_paths.sh"
STAGE="${STAGE:-A}"
R="${_REPO}/scripts/agent1_genetics/R"
PH="${AGENT1_OUT}/phenotypes"
G="${AGENT1_OUT}/gwas_summary_stats"
OUTP="${PH}/lpa_censored"
W="/path/to/scratch/cad-genetics/lpacens"
EXPORTS="${_REPO}/exports/agent1/lpa_censored"
mkdir -p "${OUTP}" "${W}/step1" "${W}/step2" "${W}/tmp" "${EXPORTS}"
RAW=/path/to/ukb_phenotypes

run_r() { set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u; Rscript "$@"; set +u; conda deactivate; set -u; }

if [[ "${STAGE}" == "A" ]]; then
  if [[ ! -s "${W}/step1/step1_pred.list" ]]; then
    echo "[$(date '+%H:%M:%S')] phenotype"
    run_r "${R}/lib_rint_censored.R" --selftest
    run_r "${R}/prep_lpa_censored.R" \
      --raw_dir "${RAW}" --keep "${PH}/regenie_keep.txt" --fam "${UKB_GENO_STEP2}.fam" \
      --covar "${PH}/covariates.txt" --pheno_q "${PH}/pheno_quantitative_statinfix.txt" \
      --lib "${R}/lib_rint_censored.R" --out "${OUTP}"

    set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_REGENIE}"; set -u
    echo "[$(date '+%H:%M:%S')] REGENIE step 1"
    regenie --step 1 \
      --bed "${UKB_GENO_STEP1}" --extract "${PH}/qc_pruned_snps.txt" \
      --keep "${PH}/regenie_keep.txt" \
      --phenoFile "${OUTP}/pheno_lpa_censored.txt" --phenoColList LPA_MEAS,LPA_CENS \
      --covarFile "${PH}/covariates.txt" \
      --bsize 1000 --lowmem --lowmem-prefix "${W}/tmp/step1" \
      --threads "${SLURM_NTASKS}" --gz --out "${W}/step1/step1"
    set +u; conda deactivate; set -u
  fi
  CHRS="6"; PAR=1; THREADS=64
else
  CHRS="1 2 3 4 5 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22"; PAR=8; THREADS=16
fi
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_REGENIE}"; set -u

step2_chr() {
  local CHR="$1" THREADS="$2" P
  P="${W}/step2/lpa_chr${CHR}"
  if [[ -s "${P}.log" ]] && grep -q "End time" "${P}.log" && [[ -s "${P}_LPA_CENS.regenie.gz" ]]; then
    echo "chr${CHR} already complete"; return 0
  fi
  regenie --step 2 \
    --bed "${UKB_GENO_STEP2}" --keep "${PH}/regenie_keep.txt" \
    --phenoFile "${OUTP}/pheno_lpa_censored.txt" --phenoColList LPA_MEAS,LPA_CENS \
    --covarFile "${PH}/covariates.txt" \
    --pred "${W}/step1/step1_pred.list" \
    --bsize 400 --chr "${CHR}" --minMAC 5 \
    --threads "${THREADS}" --gz --out "${P}" > /dev/null
  echo "chr${CHR} done $(date '+%H:%M:%S')"
}
export -f step2_chr; export UKB_GENO_STEP2 PH OUTP W
echo "[$(date '+%H:%M:%S')] REGENIE step 2: ${CHRS}"
echo ${CHRS} | tr ' ' '\n' | xargs -P "${PAR}" -I{} bash -c "step2_chr {} ${THREADS}"
set +u; conda deactivate; set -u

if [[ "${STAGE}" == "A" ]]; then
  run_r "${R}/compare_lpa_censored.R" \
    --archived "${G}/merged/lipids_LPA.txt.gz" \
    --meas "${W}/step2/lpa_chr6_LPA_MEAS.regenie.gz" --cens "${W}/step2/lpa_chr6_LPA_CENS.regenie.gz" \
    --out "${EXPORTS}/compare_chr6"
  W_TOBIT=/path/to/scratch/statint
  run_r "${R}/lpa_censored_tobit.R" \
    --input "${OUTP}/lpa_tobit_input.tsv" --dose "${W_TOBIT}/dose_LPA.raw" --weights "${W_TOBIT}/lpa_score.txt" \
    --covar "${PH}/covariates.txt" --statin "${W_TOBIT}/statin.tsv" --cores 32 --out "${EXPORTS}/tobit"
  echo STAGE_A_DONE
  exit 0
fi

# ---- Stage B: merge, then one sumstat dir per arm ----
merge_arm() {  # $1 regenie trait, $2 output file
  local T="$1" OUT="$2" tmp="$2.tmp"
  gunzip -c "${W}/step2/lpa_chr1_${T}.regenie.gz" | head -1 > "${tmp}" || true
  for c in $(seq 1 22); do gunzip -c "${W}/step2/lpa_chr${c}_${T}.regenie.gz" | tail -n +2 >> "${tmp}"; done
  gzip -f "${tmp}"; mv "${tmp}.gz" "${OUT}"
}
make_dir() {  # $1 dir: everything linked from merged_statinfix except lipids_LPA
  mkdir -p "$1"
  for f in "${G}/merged_statinfix"/*.txt.gz; do
    b=$(basename "$f"); [[ "$b" == "lipids_LPA.txt.gz" ]] && continue
    ln -sfn "$(readlink -f "$f")" "$1/$b"
  done
}
for D in merged_lpacens merged_lpacens_fixediv merged_lpameasctrl; do make_dir "${G}/${D}"; done
set +o pipefail
merge_arm LPA_CENS "${G}/merged_lpacens/lipids_LPA.txt.gz"
merge_arm LPA_MEAS "${G}/merged_lpameasctrl/lipids_LPA.txt.gz"
set -o pipefail

# Published instruments with corrected effects: selection in the MR scripts reads LOG10P only, so archived LOG10P
# with censored BETA/SE keeps exactly the published instrument sets.
cat > "${W}/fixediv.R" <<'RS'
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(TRUE)
rd <- function(f) fread(cmd = paste("gunzip -c", shQuote(f)), na.strings = "NA")
arch <- rd(a[1]); cens <- rd(a[2])
idx <- match(paste(arch$ID, arch$ALLELE0, arch$ALLELE1), paste(cens$ID, cens$ALLELE0, cens$ALLELE1))
message(sprintf("[fixediv] %d archived, %d unmatched", nrow(arch), sum(is.na(idx))))
arch[, `:=`(A1FREQ = cens$A1FREQ[idx], N = cens$N[idx], BETA = cens$BETA[idx], SE = cens$SE[idx])]
fwrite(arch, a[3], sep = " ", na = "NA", quote = FALSE, compress = "gzip")
RS
run_r "${W}/fixediv.R" "${G}/merged/lipids_LPA.txt.gz" "${G}/merged_lpacens/lipids_LPA.txt.gz" "${G}/merged_lpacens_fixediv/lipids_LPA.txt.gz"

run_r "${R}/compare_lpa_censored.R" \
  --archived "${G}/merged/lipids_LPA.txt.gz" \
  --meas "${G}/merged_lpameasctrl/lipids_LPA.txt.gz" --cens "${G}/merged_lpacens/lipids_LPA.txt.gz" \
  --out "${EXPORTS}/compare_genome"
echo STAGE_B_DONE
