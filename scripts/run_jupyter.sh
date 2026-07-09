#!/usr/bin/env bash
# ============================================================
# Jupyter Terminal Runner — Agent 1 Pipeline
# Run this instead of run_all.sh when you are already on a
# compute node (Jupyter session, idev, etc.) and do NOT want
# to submit sbatch jobs.
#
# Usage:
#   bash scripts/run_jupyter.sh              # run all steps
#   bash scripts/run_jupyter.sh --start-at 3 # resume from step
#   bash scripts/run_jupyter.sh --step 1     # run one step only
#
# Parallelism: array jobs (GWAS, coloc, MR, PRS) run in
# background batches controlled by MAX_PARALLEL.
# ============================================================

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${REPO_DIR}/config/tacc_paths.sh"

# ── tunables ────────────────────────────────────────────────
# nproc lies under SLURM single-task binding; /proc/cpuinfo is ground truth
NPROC=$(grep -c "^processor" /proc/cpuinfo 2>/dev/null || nproc --all 2>/dev/null || nproc)
THREADS_PER_JOB=6                 # threads per GWAS job (128 cores / 22 chr ≈ 5-6)
THREADS_SMALL=4                   # for coloc / MR (light)
MAX_PARALLEL=$(( NPROC / THREADS_PER_JOB ))   # GWAS jobs running at once
[[ ${MAX_PARALLEL} -lt 1 ]] && MAX_PARALLEL=1
[[ ${MAX_PARALLEL} -gt 22 ]] && MAX_PARALLEL=22  # cap at 22 (one per chromosome)

START_AT=1
ONLY_STEP=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --start-at) START_AT="$2"; shift 2 ;;
    --step)     ONLY_STEP="$2"; shift 2 ;;
    *) echo "Unknown arg: $1"; exit 1 ;;
  esac
done

SLURM_DIR="${REPO_DIR}/scripts/agent1_genetics/slurm"
LOG_DIR="${LOGS_DIR}"
mkdir -p "${LOG_DIR}"

RUN_LOG="${LOG_DIR}/jupyter_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${RUN_LOG}") 2>&1

echo "==========================================="
echo " Agent 1 — Jupyter Runner"
echo " Node: $(hostname)   Cores: ${NPROC}"
echo " MAX_PARALLEL (GWAS): ${MAX_PARALLEL}"
echo " Start at step: ${START_AT}"
echo " Started: $(date)"
echo "==========================================="

# ── helpers ─────────────────────────────────────────────────

# run_direct <script> [var=value ...]
# Runs a SLURM script as a plain bash script with given env vars.
run_direct() {
  local script="$1"; shift
  local env_prefix=""
  for kv in "$@"; do env_prefix="${env_prefix} ${kv}"; done
  echo "[$(date '+%H:%M:%S')] Running: $(basename ${script}) ${env_prefix}"
  env SLURM_NTASKS="${NPROC}" ${env_prefix} bash "${script}"
}

# run_array <script> <indices> <threads_per_job>
# Runs an array job by looping over indices, MAX_PARALLEL at a time.
run_array() {
  local script="$1"
  local indices="$2"      # e.g. "1 2 3 ... 22" or "0 1 2 ... 48"
  local tpj="${3:-${THREADS_PER_JOB}}"
  local max_par=$(( NPROC / tpj ))
  [[ ${max_par} -lt 1 ]] && max_par=1
  [[ ${max_par} -gt 16 ]] && max_par=16

  local pids=()
  local count=0

  echo "[$(date '+%H:%M:%S')] Array: $(basename ${script}) — indices: ${indices}"
  for idx in ${indices}; do
    SLURM_ARRAY_TASK_ID="${idx}" SLURM_NTASKS="${tpj}" bash "${script}" \
      > "${LOG_DIR}/array_$(basename ${script} .slurm)_${idx}_$(date +%Y%m%d).log" 2>&1 &
    pids+=($!)
    count=$(( count + 1 ))
    # throttle: wait for a slot before launching more
    if (( count % max_par == 0 )); then
      echo "[$(date '+%H:%M:%S')]   waiting for batch of ${max_par} to finish..."
      for pid in "${pids[@]}"; do wait "${pid}" || echo "  WARN: job pid ${pid} exited non-zero"; done
      pids=()
    fi
  done
  # wait for remaining
  for pid in "${pids[@]}"; do wait "${pid}" || echo "  WARN: job pid ${pid} exited non-zero"; done
  echo "[$(date '+%H:%M:%S')]   Array complete: $(basename ${script})"
}

should_run() { [[ -z "${ONLY_STEP}" && "$1" -ge "${START_AT}" ]] || [[ "${ONLY_STEP}" == "$1" ]]; }

# ── Step 1: Phenotype preparation ───────────────────────────
if should_run 1; then
  echo ""
  echo "━━━ Step 1: Phenotype preparation ━━━"
  run_direct "${SLURM_DIR}/01_prep_phenotypes.slurm"
  echo "  Outputs:"
  ls "${AGENT1_OUT}/phenotypes/" 2>/dev/null || echo "  WARNING: phenotypes dir not found"
fi

# ── Step 2: REGENIE Step 1 (null model) ─────────────────────
if should_run 2; then
  echo ""
  echo "━━━ Step 2: REGENIE Step 1 — null model (~24 h) ━━━"
  run_direct "${SLURM_DIR}/02_regenie_step1.slurm" SLURM_NTASKS="${NPROC}"
fi

# ── Steps 3 & 4: GWAS (per-chromosome arrays) ───────────────
if should_run 3; then
  echo ""
  echo "━━━ Step 3: GWAS lipids (22 chromosomes, ${MAX_PARALLEL} at a time) ━━━"
  run_array "${SLURM_DIR}/03_regenie_step2_lipids.slurm" "$(seq 1 22)" "${THREADS_PER_JOB}"
  echo ""
  echo "━━━ Step 4: GWAS outcomes (22 chromosomes) ━━━"
  run_array "${SLURM_DIR}/04_regenie_step2_outcomes.slurm" "$(seq 1 22)" "${THREADS_PER_JOB}"
fi

# ── Step 5: Merge + munge ────────────────────────────────────
if should_run 5; then
  echo ""
  echo "━━━ Step 5: Merge chromosomes + munge for LDSC ━━━"
  run_direct "${SLURM_DIR}/05_merge_and_munge.slurm"
fi

# ── Step 6: LDSC rg ─────────────────────────────────────────
if should_run 6; then
  echo ""
  echo "━━━ Step 6: LDSC genetic correlations ━━━"
  run_direct "${SLURM_DIR}/06_ldsc_rg.slurm"
fi

# ── Step 7: Colocalization (49 pairs) ───────────────────────
if should_run 7; then
  echo ""
  echo "━━━ Step 7: Colocalization — 49 pairs ━━━"
  run_array "${SLURM_DIR}/07_coloc.slurm" "$(seq 0 48)" "${THREADS_SMALL}"
fi

# ── Step 8: MR (49 pairs) ───────────────────────────────────
if should_run 8; then
  echo ""
  echo "━━━ Step 8: Mendelian randomization — 49 pairs ━━━"
  run_array "${SLURM_DIR}/08_mr.slurm" "$(seq 0 48)" "${THREADS_SMALL}"
fi

# ── Step 8b: MVMR (7 outcomes, all lipids joint) ─────────────
if should_run 8; then
  echo ""
  echo "━━━ Step 8b: Multivariable MR — 7 lipid traits jointly → 7 outcomes ━━━"
  run_array "${SLURM_DIR}/08b_mvmr.slurm" "$(seq 0 6)" "${THREADS_SMALL}"
fi

# ── Step 9: Fine-mapping ─────────────────────────────────────
if should_run 9; then
  echo ""
  echo "━━━ Step 9: SuSiE fine-mapping ━━━"
  run_direct "${SLURM_DIR}/09_finemapping.slurm"
fi

# ── Step 10: PRS (7 traits) ──────────────────────────────────
if should_run 10; then
  echo ""
  echo "━━━ Step 10: PRS — 7 lipid traits ━━━"
  run_array "${SLURM_DIR}/10_prs.slurm" "$(seq 0 6)" "${THREADS_PER_JOB}"
fi

echo ""
echo "==========================================="
echo " All steps complete: $(date)"
echo " Results: ${AGENT1_OUT}/"
echo "==========================================="
