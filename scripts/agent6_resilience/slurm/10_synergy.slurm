#!/bin/bash
#------------------------------------------------------------
# Agent 6 - synergy follow-up (preregistered S1-S3, commit 33f3252; deviation 5d2bf6a):
#   10a leads -> statin-free re-test of leads -> partition (lipid regions, split PRS-CS weights) -> PLINK2 component scores ->
#   10a assemble (L + N = total PRS check) -> 10 (S1a, S1b, S1c, S2) -> 11 (S3, full cohort).
# Submit FROM THE REPO ROOT (SENS_RUN_DIR = the statin-inclusive rerun, for S1b):
#   sbatch --export=ALL,RUN_DIR=$PWD/results/agent6_resilience/run_cad_ukbfree_v2,SENS_RUN_DIR=$PWD/results/agent6_resilience/run_cad_ukbfree_sens_statin \
#          scripts/agent6_resilience/slurm/10_synergy.slurm
#------------------------------------------------------------
#SBATCH -J a6_synergy
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o logs/a6_synergy_%j.out
#SBATCH -e logs/a6_synergy_%j.err
set -euo pipefail
_REPO="/path/to/cad-genetics"
source "${_REPO}/config/tacc_paths.sh"
: "${RUN_DIR:?export RUN_DIR (citable pipeline run directory)}"
HERE="${_REPO}/scripts/agent6_resilience"
OUT="${RUN_DIR}/10_synergy"; mkdir -p "${OUT}"
PRS_DIR="${RESULTS_DIR}/agent6_resilience/prs"
SUBSET="${TACC_SCRATCH}/agent6_prs_subset/ukb400k_prs_snps"
APOB="${AGENT1_OUT}/gwas_summary_stats/merged/lipids_APOB.txt.gz"
for f in "${APOB}" "${PRS_DIR}/cad_ukbfree_all_chr.txt" "${PRS_DIR}/prs_cad_ukbfree.txt" "${SUBSET}.bed" "${PLINK2}"; do
  [[ -e "$f" ]] || { echo "ERROR: missing $f"; exit 1; }
done
if [[ -n "${SENS_RUN_DIR:-}" ]]; then SENS=(--sens_run_dir "${SENS_RUN_DIR}" --config_sens "${HERE}/config_sens_statin.yaml"); else SENS=(); fi
CS=(--config_syn "${HERE}/config_synergy.yaml")

set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
echo "[$(date)] ApoB leads"
Rscript "${HERE}/R/10a_prs_partition.R" --step leads "${CS[@]}" --apob_gwas "${APOB}" --out "${OUT}"
echo "[$(date)] statin-free re-test of leads (deviation 2026-09-29)"
"${PLINK2}" --bfile "${UKB_GENO_STEP2}" --keep "${UKB_SAMPLE_EID}" --extract "${OUT}/apob_leads.snps" --export A \
  --threads 16 --memory 100000 --out "${OUT}/apob_leads_geno"
Rscript "${HERE}/R/10a_prs_partition.R" --step retest "${CS[@]}" --raw "${OUT}/apob_leads_geno.raw" --cohort "${RUN_DIR}/01/cohort.tsv" \
  --pheno_dir "${UKB_PHENO_DIR}" --out "${OUT}"
echo "[$(date)] partition"
Rscript "${HERE}/R/10a_prs_partition.R" --step partition "${CS[@]}" --weights "${PRS_DIR}/cad_ukbfree_all_chr.txt" --out "${OUT}"
awk 'NR > 1 {print $1}' "${OUT}/prs_component_weights.tsv" > "${OUT}/component.snps"
echo "[$(date)] score components"
# the weights file HAS a header (SNP A1 L250 N250 L500 N500 L1000 N1000), so header-read is correct here
"${PLINK2}" --bfile "${SUBSET}" --keep "${UKB_SAMPLE_EID}" --extract "${OUT}/component.snps" \
  --score "${OUT}/prs_component_weights.tsv" 1 2 header-read cols=+scoresums --score-col-nums 3-8 \
  --threads 16 --memory 100000 --out "${OUT}/prs_components"
echo "[$(date)] assemble"
Rscript "${HERE}/R/10a_prs_partition.R" --step assemble "${CS[@]}" --sscore "${OUT}/prs_components.sscore" --total_prs "${PRS_DIR}/prs_cad_ukbfree.txt" --out "${OUT}"
echo "[$(date)] stage 10"
Rscript "${HERE}/R/10_synergy.R" --config "${HERE}/config.yaml" "${CS[@]}" --run_dir "${RUN_DIR}" --components "${OUT}/prs_components.tsv" \
  --pheno_dir "${UKB_PHENO_DIR}" --out "${OUT}" "${SENS[@]}"
echo "[$(date)] stage 11"
Rscript "${HERE}/R/11_genetic_synergy.R" --config "${HERE}/config.yaml" "${CS[@]}" --run_dir "${RUN_DIR}" --components "${OUT}/prs_components.tsv" --out "${OUT}"
echo "[$(date)] done: ${OUT}"
