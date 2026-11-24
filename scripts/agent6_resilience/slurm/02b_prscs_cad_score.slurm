#!/bin/bash
#------------------------------------------------------------
# Agent 6 - merge PRS-CS posterior weights and score UKB with PLINK2.
# Output: results/agent6_resilience/prs/prs_${PREFIX}.txt  (FID IID PRS); PREFIX=cad_ukbfree (default) or cad_nikpay
# Submit FROM THE REPO ROOT after the array finishes:
#   sbatch --dependency=afterok:<ARRAY_JOB_ID> scripts/agent6_resilience/slurm/02b_prscs_cad_score.slurm
#------------------------------------------------------------
#SBATCH -J a6_prscs_score
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 06:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o logs/a6_prscs_score_%j.out
#SBATCH -e logs/a6_prscs_score_%j.err
set -euo pipefail
_REPO="/path/to/cad-genetics"
source "${_REPO}/config/tacc_paths.sh"
OUT_DIR="${RESULTS_DIR}/agent6_resilience/prs"
PREFIX="${PREFIX:-cad_ukbfree}"
W="${OUT_DIR}/${PREFIX}_all_chr.txt"
n_files=$(ls "${OUT_DIR}"/${PREFIX}_pst_eff_a1_b0.5_phiauto_chr*.txt | wc -l)
[[ "${n_files}" -eq 22 ]] || { echo "ERROR: expected 22 per-chr weight files, found ${n_files}"; exit 1; }
cat "${OUT_DIR}"/${PREFIX}_pst_eff_a1_b0.5_phiauto_chr{1..22}.txt > "${W}"
echo "weights: $(wc -l < "${W}") SNPs"
# Restrict PLINK2 to the weighted SNPs: without --extract it scans all 14.8M
# variants x 402k samples for allele frequencies and exceeds the walltime.
awk '{print $2}' "${W}" > "${W}.snps"
# PRS-CS columns: CHR SNP BP A1 A2 BETA  -> score on SNP(2) A1(4) BETA(6).
# PRS-CS output has NO header line: never pass header-read (it silently eats the first SNP
# and renames the score column after its beta).
# Prefer the scratch subset bed built by 02c_extract_prs_bed.slurm (minutes instead of hours)
SUBSET="${TACC_SCRATCH}/agent6_prs_subset/ukb400k_prs_snps"
if [[ -s "${SUBSET}.bed" ]]; then BFILE="${SUBSET}"; echo "using subset bed ${SUBSET}"; else BFILE="${UKB_GENO_STEP2}"; echo "WARNING: subset bed absent, streaming full corral bed"; fi
"${PLINK2}" --bfile "${BFILE}" --keep "${UKB_SAMPLE_EID}" --extract "${W}.snps" \
  --score "${W}" 2 4 6 cols=+scoresums list-variants --threads 16 --memory 200000 \
  --out "${OUT_DIR}/${PREFIX}_score"
n_used=$(wc -l < "${OUT_DIR}/${PREFIX}_score.sscore.vars")
echo "variants used in score: ${n_used} of $(wc -l < "${W}.snps")"
[[ "${n_used}" -gt 500000 ]] || { echo "ERROR: fewer than 500k weighted SNPs matched the genotype file"; exit 1; }
awk 'BEGIN{OFS="\t"} NR==1{for(i=1;i<=NF;i++){if($i=="#FID"||$i=="FID")f=i;if($i=="IID")id=i;if($i=="SCORE1_SUM")s=i}; if(!f||!id||!s){print "ERROR: FID/IID/SCORE1_SUM column missing in sscore header: " $0 > "/dev/stderr"; exit 1}; print "FID","IID","PRS"; next} {print $f,$id,$s}' \
  "${OUT_DIR}/${PREFIX}_score.sscore" > "${OUT_DIR}/prs_${PREFIX}.txt"
bad=$(awk -F'\t' 'NF!=3' "${OUT_DIR}/prs_${PREFIX}.txt" | wc -l)
[[ "${bad}" -eq 0 ]] || { echo "ERROR: ${bad} rows in PRS file do not have exactly 3 columns"; exit 1; }
echo "scored: $(($(wc -l < "${OUT_DIR}/prs_${PREFIX}.txt") - 1)) individuals"
md5sum "${W}" "${OUT_DIR}/prs_${PREFIX}.txt" > "${OUT_DIR}/prs_${PREFIX}.md5"
