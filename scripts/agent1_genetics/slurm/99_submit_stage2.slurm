#!/bin/bash
#------------------------------------------------------------
# Stage 2 launcher: submits steps 6-10 after step 5 (merge/munge) completes.
#
# Deferred submission avoids QOSMaxSubmitJobPerUserLimit on TACC:
# steps 3+4 GWAS arrays (44 tasks) and mask/burden arrays (44+ tasks)
# are all pending at pipeline launch. By the time step 5 finishes,
# those arrays are done and quota is free for the large coloc/MR arrays.
#
# Submitted automatically by run_all.sh with --dependency=afterok:JOB5.
# RUN_LOG_DIR and RUN_ID are inherited via --export=ALL from run_all.sh.
#------------------------------------------------------------
#SBATCH -J cad_stage2_launcher
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:15:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=fail
#SBATCH --mail-user=YOUR_EMAIL

set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

# RUN_LOG_DIR exported by run_all.sh; fall back if launched standalone
RUN_LOG_DIR="${RUN_LOG_DIR:-${LOGS_DIR}/runs/stage2_$(date +%Y%m%d_%H%M%S)}"
mkdir -p "${RUN_LOG_DIR}"

SLURM_DIR="${_REPO}/scripts/agent1_genetics/slurm"

echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME}/${SLURM_JOB_ID}] Stage 2 launcher — submitting steps 6-10"
echo "  Log dir: ${RUN_LOG_DIR}"

# Inline submit helper (mirrors run_all.sh submit())
_submit() {
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
    --chdir="${SLURM_DIR}" \
    "${SLURM_DIR}/${script}" 2>&1)
  jid=$(echo "${raw}" | grep -oP '^\d+$' | tail -1)
  if [[ -z "${jid}" ]]; then
    echo "  ERROR: sbatch failed for ${script}" >&2
    echo "  sbatch output: ${raw}" >&2
    exit 1
  fi
  echo "  Submitted ${script}: job ${jid}  [→ ${step_name}]"
  printf '%s\t%s\t%s\t%s\n' "${step_name}" "${script}" "${jid}" "$(date +%Y-%m-%dT%H:%M:%S)" \
    >> "${RUN_LOG_DIR}/manifest.tsv"
  echo "${jid}"
}

# Steps 6, 7, 8, 8b, 10 all run in parallel (no inter-dependencies).
# Step 9 waits for coloc (step 7).
JOB6=$(_submit  06_ldsc_rg.slurm)
JOB7=$(_submit  07_coloc.slurm      --array "0-62")
JOB8=$(_submit  08_mr.slurm         --array "0-62")
JOB8B=$(_submit 08b_mvmr.slurm      --array "0-6")
JOB9=$(_submit  09_finemapping.slurm --dependency "${JOB7}")
JOB10=$(_submit 10_prs.slurm        --array "0-6")

echo ""
echo "  Step 6  LDSC rg:      ${JOB6}"
echo "  Step 7  Coloc ×63:    ${JOB7}"
echo "  Step 8  MR ×63:       ${JOB8}"
echo "  Step 8b MVMR ×7:      ${JOB8B}"
echo "  Step 9  Finemapping:  ${JOB9}  (after coloc ${JOB7})"
echo "  Step 10 PRS ×7:       ${JOB10}"
echo ""
echo "[$(date '+%H:%M:%S')] Stage 2 launcher complete."
