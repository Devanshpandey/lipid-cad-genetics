#!/bin/bash
#------------------------------------------------------------
# SLURM: merge per-chromosome burden annotation files
# Run after prepare_burden_masks.slurm (array 1-22) completes.
#
# Reads:  ${UKB_BURDEN_ANNO_DIR}/by_chr/anno_chr{1-22}.txt
#         ${UKB_BURDEN_ANNO_DIR}/by_chr/setlist_chr{1-22}.txt
# Writes: ${UKB_BURDEN_ANNO_FILE}   (REGENIE --anno-file)
#         ${UKB_BURDEN_SETLIST}     (REGENIE --set-list)
#         ${UKB_BURDEN_MASKS}       (REGENIE --mask-def)
#
# Submit: sbatch merge_burden_masks.slurm
#   (or: sbatch --dependency=afterok:<array_jobid> merge_burden_masks.slurm)
#------------------------------------------------------------
#SBATCH -J burden_merge
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH -A OTH26002
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=devanshpandey@utexas.edu

set -euo pipefail
# Resolve repo root: SLURM_SUBMIT_DIR when run via sbatch, hardcoded path as fallback
# (BASH_SOURCE[0] points to /var/spool in SLURM jobs — cannot use dirname approach)
_REPO="${SLURM_SUBMIT_DIR:-/work/07880/devansh/lonestar/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_PY}"
set -u

BY_CHR_DIR="${UKB_BURDEN_ANNO_DIR}/by_chr"
echo "[$(date '+%H:%M:%S')] [${SLURM_JOB_NAME:-shell}/${SLURM_JOB_ID:-0}] Merging per-chromosome burden annotation files..."

# ---- Check all 22 chromosomes are present ----
MISSING=0
for CHR in $(seq 1 22); do
  f="${BY_CHR_DIR}/anno_chr${CHR}.txt"
  if [[ ! -f "${f}" || ! -s "${f}" ]]; then
    echo "  WARNING: missing or empty: ${f}"
    MISSING=$((MISSING + 1))
  fi
done
if [[ "${MISSING}" -gt 0 ]]; then
  echo "ERROR: ${MISSING} chromosome(s) missing — check prepare_burden_masks.slurm logs"
  exit 1
fi
echo "  All 22 chromosomes present."

# ---- Merge anno files (cat, no header) ----
echo "[$(date)] Merging anno files..."
: > "${UKB_BURDEN_ANNO_FILE}"
for CHR in $(seq 1 22); do
  cat "${BY_CHR_DIR}/anno_chr${CHR}.txt" >> "${UKB_BURDEN_ANNO_FILE}"
done
N_ANNO=$(wc -l < "${UKB_BURDEN_ANNO_FILE}")
echo "  Total annotated variants: ${N_ANNO}"

# ---- Merge set list files ----
echo "[$(date)] Merging set list files..."
: > "${UKB_BURDEN_SETLIST}"
for CHR in $(seq 1 22); do
  f="${BY_CHR_DIR}/setlist_chr${CHR}.txt"
  [[ -f "${f}" ]] && cat "${f}" >> "${UKB_BURDEN_SETLIST}"
done
N_GENES=$(wc -l < "${UKB_BURDEN_SETLIST}")
echo "  Total genes: ${N_GENES}"

# ---- Write mask definition file ----
cat > "${UKB_BURDEN_MASKS}" << 'EOF'
mask1	pLoF
mask2	pLoF,missense3
mask3	pLoF,missense3,missense2
EOF
echo "[$(date)] Mask definitions:"
cat "${UKB_BURDEN_MASKS}"

# ---- Annotation category summary ----
echo ""
echo "[$(date)] === Annotation category counts ==="
awk '{print $3}' "${UKB_BURDEN_ANNO_FILE}" | sort | uniq -c | sort -rn
echo ""

# ---- Validation ----
echo "[$(date)] === File summary ==="
for label in "anno_file:${UKB_BURDEN_ANNO_FILE}" \
             "set_list:${UKB_BURDEN_SETLIST}" \
             "mask_def:${UKB_BURDEN_MASKS}"; do
  name="${label%%:*}"
  path="${label#*:}"
  if [[ -f "${path}" && -s "${path}" ]]; then
    echo "  OK  ${name}: $(wc -l < ${path}) lines  →  ${path}"
  else
    echo "  FAIL  ${name}: ${path}"
  fi
done

echo ""
echo "[$(date)] Burden mask files ready."
echo "Next step: sbatch --array=1-22 scripts/agent1_genetics/slurm/03b_burden_test.slurm"
