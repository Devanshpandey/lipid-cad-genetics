#!/usr/bin/env bash
# ============================================================
# Agent 1 — Genetics Pipeline Orchestrator  (v2 — 2026-06)
# Run on a Lonestar6 LOGIN NODE.
# Submits SLURM jobs in dependency order.
#
# Changes vs v1:
#   Step 3b — REGENIE gene-based burden testing (exome, additive+recessive)
#   Steps 7/8/LDSC — 63 pairs (7 lipids × 9 outcomes; added AF, PAD)
#   Phenotypes — nonHDL_C, eGFR, CRP, HBA1C; fixed MACE; REVASC+OPCS; STATIN_USE
#
# Usage:
#   bash scripts/agent1_genetics/run_all.sh
#   bash scripts/agent1_genetics/run_all.sh --start-at 3   # resume from step 3
#   bash scripts/agent1_genetics/run_all.sh --skip-burden  # skip exome burden (if no WES data)
# ============================================================

set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${REPO_DIR}/config/tacc_paths.sh"

START_AT=1
SKIP_BURDEN=0
BUILD_MASKS=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --start-at)    START_AT="$2"; shift 2 ;;
    --skip-burden) SKIP_BURDEN=1; shift ;;
    --build-masks) BUILD_MASKS=1; shift ;;   # submit mask prep jobs before step 3b
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

SLURM_DIR="${REPO_DIR}/scripts/agent1_genetics/slurm"
RUN_ID="$(date +%Y%m%d_%H%M%S)"
RUN_LOG_DIR="${LOGS_DIR}/runs/${RUN_ID}"
mkdir -p "${RUN_LOG_DIR}"
export RUN_ID RUN_LOG_DIR
LOG_FILE="${RUN_LOG_DIR}/run_all.log"
exec > >(tee -a "${LOG_FILE}") 2>&1
printf 'step\tscript\tjob_id\tsubmitted_at\n' > "${RUN_LOG_DIR}/manifest.tsv"

echo "================================================"
echo " Agent 1 — Genetics Pipeline  (v2)"
echo " Started: $(date)"
echo " TACC_WORK: ${TACC_WORK}"
echo " Starting at step: ${START_AT}"
echo " Skip burden test: ${SKIP_BURDEN}"
echo "================================================"

# ---- submit() helper ----
# submit <slurm_script> [--dependency <jobid>] [--array <spec>]
# Returns only the numeric SLURM job ID on stdout.
# All human-readable messages go to stderr so they don't pollute $() captures.
# Uses grep to strip TACC's MOTD banner from sbatch --parsable output.
submit() {
  local script="$1"; shift
  local dep="" arr="" step_name="" is_array=0
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --dependency) dep="--dependency=afterok:$2"; shift 2 ;;
      --array)      arr="--array=$2"; is_array=1; shift 2 ;;
      --name)       step_name="$2"; shift 2 ;;
      *) shift ;;
    esac
  done
  step_name="${step_name:-${script%.slurm}}"

  local out_pat err_pat
  if [[ "${is_array}" -eq 1 ]]; then
    out_pat="${RUN_LOG_DIR}/${step_name}_%A_%a.out"
    err_pat="${RUN_LOG_DIR}/${step_name}_%A_%a.err"
  else
    out_pat="${RUN_LOG_DIR}/${step_name}_%j.out"
    err_pat="${RUN_LOG_DIR}/${step_name}_%j.err"
  fi

  local raw jid
  raw=$(sbatch --parsable ${dep} ${arr} \
    --output="${out_pat}" --error="${err_pat}" \
    --chdir="${REPO_DIR}/scripts/agent1_genetics/slurm" \
    "${SLURM_DIR}/${script}" 2>&1)
  # Extract purely numeric job ID — strips TACC MOTD banner lines
  jid=$(echo "${raw}" | grep -oP '^\d+$' | tail -1)
  if [[ -z "${jid}" ]]; then
    echo "  ERROR: sbatch returned no job ID for ${script}" >&2
    echo "  sbatch output: ${raw}" >&2
    exit 1
  fi
  echo "  Submitted ${script}: job ${jid}  [→ ${step_name}]" >&2
  printf '%s\t%s\t%s\t%s\n' "${step_name}" "${script}" "${jid}" "$(date +%Y-%m-%dT%H:%M:%S)" \
    >> "${RUN_LOG_DIR}/manifest.tsv"
  echo "${jid}"
}

# ---- multi-dependency helper ----
# join_deps "j1 j2 j3" → "j1:j2:j3"
join_deps() {
  echo "$1" | tr ' ' ':'
}

check_inputs() {
  echo "[check] Verifying required inputs..."
  [[ -d "${UKB_PHENO_DIR}" ]]       || { echo "ERROR: UKB_PHENO_DIR not found: ${UKB_PHENO_DIR}"; exit 1; }
  [[ -f "${UKB_SAMPLE_EID}" ]]      || { echo "ERROR: UKB_SAMPLE_EID not found: ${UKB_SAMPLE_EID}"; exit 1; }
  [[ -f "${UKB_GENO_STEP1}.bed" ]]  || { echo "ERROR: GENO_STEP1 BED not found: ${UKB_GENO_STEP1}.bed"; exit 1; }
  [[ -f "${UKB_GENO_STEP2}.bed" ]]  || { echo "ERROR: GENO_STEP2 BED not found: ${UKB_GENO_STEP2}.bed"; exit 1; }
  if [[ "${SKIP_BURDEN}" -eq 0 ]]; then
    if [[ ! -f "${UKB_BURDEN_ANNO_FILE}" ]]; then
      if [[ "${BUILD_MASKS}" -eq 1 ]]; then
        echo "[check] Burden annotation file absent — mask prep jobs will be submitted automatically (--build-masks)."
      else
        echo ""
        echo "WARNING: Burden annotation file not found: ${UKB_BURDEN_ANNO_FILE}"
        echo "  Options:"
        echo "    --build-masks  : auto-submit mask prep jobs before the burden test (recommended for fresh run)"
        echo "    --skip-burden  : skip step 3b entirely"
        echo ""
        read -r -p "  Auto-submit mask prep and continue with burden test? (y/N): " reply
        if [[ "${reply,,}" == "y" ]]; then
          BUILD_MASKS=1
        else
          echo "  Re-run with --build-masks or --skip-burden."
          exit 1
        fi
      fi
    fi
  fi
  echo "[check] All inputs present."
}

# ============================================================
# Step 1: Phenotype preparation
# ============================================================
if [[ "${START_AT}" -le 1 ]]; then
  check_inputs
  echo "[step 1] Phenotype preparation..."
  JOB1=$(submit 01_prep_phenotypes.slurm)
else
  JOB1=""
  echo "[skip] Step 1 — assuming phenotypes already prepared"
fi

# ============================================================
# Step 2: REGENIE Step 1 (null model — common variants)
# ============================================================
if [[ "${START_AT}" -le 2 ]]; then
  echo "[step 2] REGENIE Step 1 (null model)..."
  if [[ -n "${JOB1}" ]]; then
    JOB2=$(submit 02_regenie_step1.slurm --dependency "${JOB1}")
  else
    JOB2=$(submit 02_regenie_step1.slurm)
  fi
else
  JOB2=""
  echo "[skip] Step 2 — assuming REGENIE Step 1 already done"
fi

# ============================================================
# Steps 3 + 3b + 4: GWAS (lipids chr1-22), burden test (chr1-22),
#                    GWAS (outcomes chr1-22) — all in parallel
# ============================================================
if [[ "${START_AT}" -le 3 ]]; then
  echo "[step 3] REGENIE Step 2 — Lipid GWAS (chr1-22)..."
  if [[ -n "${JOB2}" ]]; then
    JOB3=$(submit 03_regenie_step2_lipids.slurm   --dependency "${JOB2}" --array "1-22")
    JOB4=$(submit 04_regenie_step2_outcomes.slurm --dependency "${JOB2}" --array "1-22")
  else
    JOB3=$(submit 03_regenie_step2_lipids.slurm   --array "1-22")
    JOB4=$(submit 04_regenie_step2_outcomes.slurm --array "1-22")
  fi
  echo "[step 4] REGENIE Step 2 — Outcome GWAS (chr1-22)... job ${JOB4}"

  if [[ "${SKIP_BURDEN}" -eq 0 ]]; then
    # Optionally build annotation masks first (if not already built)
    MASK_DEP=""
    if [[ "${BUILD_MASKS}" -eq 1 || ! -f "${UKB_BURDEN_ANNO_FILE}" ]]; then
      echo "[step 3b-prep] Submitting burden mask preparation (chr1-22 array + merge)..."
      JOB_MASK_ARRAY=$(submit prepare_burden_masks.slurm --array "1-22")
      JOB_MASK_MERGE=$(submit merge_burden_masks.slurm --dependency "${JOB_MASK_ARRAY}")
      MASK_DEP="${JOB_MASK_MERGE}"
      echo "  Mask prep array: ${JOB_MASK_ARRAY}"
      echo "  Mask merge job:  ${JOB_MASK_MERGE}"
    else
      echo "[step 3b-prep] Annotation files already exist — skipping mask prep"
    fi

    echo "[step 3b] REGENIE burden test — WES (chr1-22, additive+recessive)..."
    BURDEN_DEP=""
    [[ -n "${JOB2}" ]]      && BURDEN_DEP="${JOB2}"
    [[ -n "${MASK_DEP}" ]]  && BURDEN_DEP="${BURDEN_DEP:+${BURDEN_DEP}:}${MASK_DEP}"

    if [[ -n "${BURDEN_DEP}" ]]; then
      JOB3B=$(submit 03b_burden_test.slurm --dependency "${BURDEN_DEP}" --array "1-22")
    else
      JOB3B=$(submit 03b_burden_test.slurm --array "1-22")
    fi
    echo "  Burden test job: ${JOB3B}"
  else
    JOB3B=""
    echo "[skip] Step 3b — burden test skipped (--skip-burden)"
  fi
else
  JOB3=""; JOB4=""; JOB3B=""
  echo "[skip] Steps 3/3b/4 — assuming GWAS + burden already done"
fi

# ============================================================
# Step 5: Merge chromosomes + munge for LDSC
#         Waits for both lipid and outcome GWAS to complete
# ============================================================
if [[ "${START_AT}" -le 5 ]]; then
  echo "[step 5] Merge & munge sumstats..."
  MERGE_DEPS="${JOB3:+${JOB3}} ${JOB4:+${JOB4}}"
  MERGE_DEPS="${MERGE_DEPS// /}"  # remove spaces if either is empty
  if [[ -n "${JOB3}" && -n "${JOB4}" ]]; then
    JOB5=$(submit 05_merge_and_munge.slurm --dependency "$(join_deps "${JOB3} ${JOB4}")")
  elif [[ -n "${JOB3}" ]]; then
    JOB5=$(submit 05_merge_and_munge.slurm --dependency "${JOB3}")
  else
    JOB5=$(submit 05_merge_and_munge.slurm)
  fi
else
  JOB5=""
  echo "[skip] Step 5 — assuming merge/munge already done"
fi

# ============================================================
# Steps 6-10: LDSC, coloc, MR, MVMR, finemapping, PRS
#
# These are NOT submitted now — they go via a stage2 launcher job
# that fires after step 5 completes.  By that point the 88 GWAS
# array tasks (steps 3+4) are guaranteed done, freeing up quota for
# the large coloc/MR arrays (63+63 tasks) without hitting
# QOSMaxSubmitJobPerUserLimit.
#
# If resuming from step 6+ (--start-at 6), submit directly instead.
# ============================================================

if [[ "${START_AT}" -ge 6 ]]; then
  # Manual resume: step 5 already done, submit stage2 steps directly now.
  echo "[steps 6-10] Resuming from step ${START_AT} — submitting directly..."
  JOB6=""
  JOB7=""
  JOB8=""
  JOB8B=""
  JOB9=""

  [[ "${START_AT}" -le 6 ]] && JOB6=$(submit 06_ldsc_rg.slurm)
  [[ "${START_AT}" -le 7 ]] && JOB7=$(submit 07_coloc.slurm   --array "0-62") \
                            || echo "[skip] Step 7 — assuming coloc done"
  [[ "${START_AT}" -le 8 ]] && JOB8=$(submit 08_mr.slurm      --array "0-62") \
                            || echo "[skip] Step 8 — assuming MR done"
  [[ "${START_AT}" -le 8 ]] && JOB8B=$(submit 08b_mvmr.slurm  --array "0-6") \
                            || true

  if [[ "${START_AT}" -le 9 ]]; then
    if [[ -n "${JOB7}" ]]; then
      JOB9=$(submit 09_finemapping.slurm --dependency "${JOB7}")
    else
      JOB9=$(submit 09_finemapping.slurm)
    fi
  else
    echo "[skip] Step 9 — assuming fine-mapping done"
  fi

  [[ "${START_AT}" -le 10 ]] && submit 10_prs.slurm --array "0-6" \
                              || echo "[skip] Step 10 — assuming PRS done"

elif [[ "${START_AT}" -le 5 ]]; then
  # Normal flow: defer steps 6-10 to a launcher job after step 5.
  echo "[steps 6-10] Submitting stage2 launcher (fires after step 5 completes)..."
  if [[ -n "${JOB5}" ]]; then
    JOB_STAGE2=$(submit 99_submit_stage2.slurm --dependency "${JOB5}" --name "stage2_launcher")
  else
    JOB_STAGE2=$(submit 99_submit_stage2.slurm --name "stage2_launcher")
  fi
  echo "  Stage 2 launcher: job ${JOB_STAGE2}"
  echo "  (Steps 6-10 will be queued automatically when step 5 finishes)"
fi

# ============================================================
# Summary
# ============================================================
echo ""
echo "================================================"
echo " All Agent 1 jobs submitted."
echo " Run ID:        ${RUN_ID}"
echo " Log directory: ${RUN_LOG_DIR}/"
echo " Manifest:      ${RUN_LOG_DIR}/manifest.tsv"
echo " Monitor with:  squeue -u ${USER}"
echo " Results will appear in: ${AGENT1_OUT}/"
echo " Pull to Mac:  bash scripts/sync_results.sh agent1"
echo "================================================"

echo ""
echo "Job dependency graph (v2):"
echo ""
echo "  [01 pheno] → [02 regenie_step1]"
echo "                    ├─→ [03 gwas_lipids       ×22] ─┐"
echo "                    ├─→ [04 gwas_outcomes     ×22] ─┤→ [05 merge/munge]"
echo "                    └─→ [03b burden_test      ×22] ─┘"
echo "                         ↑ depends on:"
echo "                    [mask_prep ×22] → [mask_merge]"
echo "                    (auto-submitted if annotation files missing)"
echo ""
echo "  [05 merge] → [stage2 launcher] → [06 ldsc_rg]"
echo "                                → [07 coloc  ×63] → [09 finemapping]"
echo "                                → [08 MR     ×63]"
echo "                                → [08b MVMR  ×7]"
echo "                                → [10 PRS    ×7]"
echo "  (Stage 2 is deferred to avoid QOSMaxSubmitJobPerUserLimit"
echo "   while steps 3+4 GWAS arrays are still pending)"
echo ""
echo "WES annotation source:"
echo "  ${UKB_WES_JOYCE_VCF}/annotated_chr_{1-22}_sites_dbNSFP.vcf.gz"
if [[ "${SKIP_BURDEN}" -eq 1 ]]; then
  echo ""
  echo "  NOTE: Step 3b (burden test) was SKIPPED (--skip-burden)."
  echo "  To run later:"
  echo "    sbatch --array=1-22 scripts/agent1_genetics/slurm/prepare_burden_masks.slurm"
  echo "    # wait for completion, then:"
  echo "    sbatch scripts/agent1_genetics/slurm/merge_burden_masks.slurm"
  echo "    sbatch --array=1-22 scripts/agent1_genetics/slurm/03b_burden_test.slurm"
fi
echo ""
echo "New in v2 (vs previous run):"
echo "  • Step 3b: WES burden test (additive+recessive, mask1/mask2, pLoF+damaging missense)"
echo "  • Mask source: joyce_annotated_vcf/ dbNSFP VCFs (VEP+LOFTEE+REVEL+CADD)"
echo "  • New outcomes: AF (I48), PAD (I70/I73/I74) — 9 outcomes total"
echo "  • New quant traits: nonHDL_C, eGFR (CKD-EPI 2021), CRP, HbA1c"
echo "  • Fixed MACE: now includes CV_DEATH"
echo "  • REVASC: supplemented with OPCS4 PCI/CABG codes (fid41200/41201)"
echo "  • STATIN_USE as a binary phenotype"
