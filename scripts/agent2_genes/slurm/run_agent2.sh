#!/bin/bash
# ============================================================
# Agent 2 — Gene Prioritization
# Run via nohup on Jupyter compute node (128 cores)
#
# Usage:
#   nohup bash scripts/agent2_genes/slurm/run_agent2.sh \
#     > logs/agent2_$(date +%Y%m%d_%H%M%S).log 2>&1 &
# ============================================================
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
source "${REPO_DIR}/config/tacc_paths.sh"

LOG_FILE="${LOGS_DIR}/agent2_$(date +%Y%m%d).log"
mkdir -p "${LOGS_DIR}" "${AGENT2_OUT}"
exec > >(tee -a "${LOG_FILE}") 2>&1

# Reference data paths
REF="${TACC_WORK}/data/reference"
GENE_BED="${REF}/ensembl_genes_hg38.bed"
EQTL_DIR="${REF}/eqtl_catalogue"
OT_L2G_DIR="${REF}/opentargets/l2g"
OT_TARGET_DIR="${REF}/opentargets/target"
OT_DRUG_DIR="${REF}/opentargets/known_drug"

# Agent 1 inputs
CREDIBLE_SETS="${AGENT1_OUT}/finemapping_credible_sets.csv"
COLOC_RESULTS="${AGENT1_OUT}/coloc_results.csv"
SUMSTAT_DIR="${AGENT1_OUT}/gwas_summary_stats/merged"
MR_RESULTS="${AGENT1_OUT}/mr_all_pairs.csv"

# Agent 2 outputs
POS_OUT="${AGENT2_OUT}/positional_mapping.csv"
EQTL_OUT="${AGENT2_OUT}/gtex_eqtl_coloc.csv"
L2G_OUT="${AGENT2_OUT}/opentargets_l2g.csv"
RANKED_OUT="${AGENT2_OUT}/ranked_causal_gene_table.csv"

PYTHON_SCRIPTS="${SCRIPT_DIR}/../python"
R_SCRIPTS="${SCRIPT_DIR}/../R"

set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
set -u

echo "=============================================="
echo " Agent 2: Gene Prioritization"
echo " Started: $(date)"
echo "=============================================="

# ── Check Agent 1 inputs ──────────────────────────────────────
echo ""
echo "[$(date)] Checking Agent 1 inputs..."
for f in "${CREDIBLE_SETS}" "${COLOC_RESULTS}" "${MR_RESULTS}"; do
  if [[ ! -f "${f}" ]]; then
    echo "ERROR: Required Agent 1 output missing: ${f}"
    exit 1
  fi
done
echo "  All Agent 1 inputs present"

# ── Check reference data ──────────────────────────────────────
echo ""
echo "[$(date)] Checking reference data..."
echo "  Gene BED:    $([ -f "${GENE_BED}" ] && echo "FOUND" || echo "MISSING - will auto-download")"
echo "  GTEx eQTL:   $([ -d "${EQTL_DIR}" ] && ls "${EQTL_DIR}"/*.tsv.gz 2>/dev/null | wc -l || echo 0) tissue files"
echo "  OT L2G:      $([ -d "${OT_L2G_DIR}" ] && ls "${OT_L2G_DIR}"/*.parquet 2>/dev/null | wc -l || echo 0) parquet files"
echo "  OT Target:   $([ -d "${OT_TARGET_DIR}" ] && echo "FOUND" || echo "MISSING")"
echo "  OT Drug:     $([ -d "${OT_DRUG_DIR}" ] && echo "FOUND" || echo "MISSING")"

EQTL_AVAILABLE=false
OT_AVAILABLE=false
[[ $(ls "${EQTL_DIR}"/*.tsv.gz 2>/dev/null | wc -l) -gt 0 ]] && EQTL_AVAILABLE=true
[[ $(ls "${OT_L2G_DIR}"/*.parquet 2>/dev/null | wc -l) -gt 0 ]] && OT_AVAILABLE=true

if [[ "${EQTL_AVAILABLE}" == "false" ]] || [[ "${OT_AVAILABLE}" == "false" ]]; then
  echo ""
  echo "  WARNING: Some reference data missing."
  echo "  Run first: nohup bash ${REPO_DIR}/scripts/download_agent2_refs.sh &"
  echo "  Pipeline will continue with available data only."
fi

# ── Step 1: Positional mapping ────────────────────────────────
echo ""
echo "[$(date)] === Step 1: Positional gene mapping ==="
conda activate "${CONDA_ENV_PY}"

if [[ ! -f "${POS_OUT}" ]]; then
  python "${PYTHON_SCRIPTS}/01_positional_mapping.py" \
    --credible_sets "${CREDIBLE_SETS}" \
    --coloc_results  "${COLOC_RESULTS}" \
    --gene_bed       "${GENE_BED}" \
    --window_kb      500 \
    --out            "${POS_OUT}"
else
  echo "  Positional mapping already done — skipping"
fi
echo "[$(date)] Step 1 complete"

# ── Step 2: GTEx eQTL colocalization ─────────────────────────
echo ""
echo "[$(date)] === Step 2: GTEx eQTL colocalization ==="
conda activate "${CONDA_ENV_R}"

if [[ ! -f "${EQTL_OUT}" ]]; then
  if [[ "${EQTL_AVAILABLE}" == "true" ]]; then
    Rscript "${R_SCRIPTS}/02_gtex_eqtl_coloc.R" \
      --credible_sets "${CREDIBLE_SETS}" \
      --sumstat_dir   "${SUMSTAT_DIR}" \
      --eqtl_dir      "${EQTL_DIR}" \
      --window_kb     500 \
      --pp4_thresh    0.5 \
      --out           "${EQTL_OUT}"
  else
    echo "  GTEx eQTL data not available — creating empty placeholder"
    echo "exposure,locus_chr,locus_bp,sentinel,tissue,gene_id,gene_name,n_snps,PP.H0,PP.H1,PP.H2,PP.H3,PP.H4,eqtl_coloc_sig" > "${EQTL_OUT}"
  fi
else
  echo "  GTEx eQTL coloc already done — skipping"
fi
echo "[$(date)] Step 2 complete"

# ── Step 3: OpenTargets L2G ───────────────────────────────────
echo ""
echo "[$(date)] === Step 3: OpenTargets L2G + druggability ==="
conda activate "${CONDA_ENV_PY}"

if [[ ! -f "${L2G_OUT}" ]]; then
  python "${PYTHON_SCRIPTS}/03_opentargets_l2g.py" \
    --credible_sets  "${CREDIBLE_SETS}" \
    --ot_l2g_dir     "${OT_L2G_DIR}" \
    --ot_target_dir  "${OT_TARGET_DIR}" \
    --ot_drug_dir    "${OT_DRUG_DIR}" \
    --out            "${L2G_OUT}"
else
  echo "  OpenTargets L2G already done — skipping"
fi
echo "[$(date)] Step 3 complete"

# ── Step 4: Integrate and score ───────────────────────────────
echo ""
echo "[$(date)] === Step 4: Evidence integration and gene scoring ==="
conda activate "${CONDA_ENV_PY}"

python "${PYTHON_SCRIPTS}/04_score_genes.py" \
  --positional    "${POS_OUT}" \
  --eqtl_coloc    "${EQTL_OUT}" \
  --ot_l2g        "${L2G_OUT}" \
  --credible_sets "${CREDIBLE_SETS}" \
  --mr_results    "${MR_RESULTS}" \
  --out           "${RANKED_OUT}"

echo "[$(date)] Step 4 complete"

# ── Export key outputs ────────────────────────────────────────
echo ""
echo "[$(date)] === Exporting results ==="
EXPORT_DIR="${REPO_DIR}/exports/agent2"
mkdir -p "${EXPORT_DIR}"
for f in "${RANKED_OUT}" "${EQTL_OUT}" "${POS_OUT}" "${L2G_OUT}"; do
  [[ -f "${f}" ]] && cp "${f}" "${EXPORT_DIR}/" && echo "  Exported: $(basename ${f})"
done

echo ""
echo "=============================================="
echo " Agent 2 complete: $(date)"
echo " Outputs: ${AGENT2_OUT}/"
ls -lh "${AGENT2_OUT}/" 2>/dev/null
echo "=============================================="
