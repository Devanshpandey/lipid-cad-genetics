#!/usr/bin/env bash
# Consolidate key outputs from all agents into results/final/
# Run on TACC after all agents complete.

set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${REPO_DIR}/config/tacc_paths.sh"

FINAL="${RESULTS_DIR}/final"
mkdir -p "${FINAL}"

echo "[consolidate] Copying final outputs to ${FINAL}/"

copy_if_exists() {
  local src="$1" dst="$2"
  if [[ -f "${src}" ]]; then
    cp "${src}" "${dst}"
    echo "  ✓ $(basename ${dst})"
  else
    echo "  ✗ MISSING: ${src}"
  fi
}

# Agent 1
copy_if_exists "${AGENT1_OUT}/ldsc_rg_matrix.csv"              "${FINAL}/rg_matrix.csv"
copy_if_exists "${AGENT1_OUT}/coloc_results.csv"               "${FINAL}/coloc_results.csv"
copy_if_exists "${AGENT1_OUT}/finemapping_credible_sets.csv"   "${FINAL}/finemapping_credible_sets.csv"

# Combine all MR results into one table
MR_DIR="${AGENT1_OUT}/mr"
if [[ -d "${MR_DIR}" ]]; then
  python - <<PYEOF
import os, polars as pl
mr_dir = "${MR_DIR}"
files = [os.path.join(mr_dir, f) for f in os.listdir(mr_dir) if f.endswith("_mr.csv")]
dfs = [pl.read_csv(f) for f in files if os.path.getsize(f) > 0]
if dfs:
    combined = pl.concat(dfs)
    combined.write_csv("${FINAL}/mr_causal_evidence_table.csv")
    print(f"  ✓ mr_causal_evidence_table.csv ({len(combined)} rows)")
PYEOF
fi

# Agent 2
copy_if_exists "${AGENT2_OUT}/ranked_causal_gene_table.csv"  "${FINAL}/ranked_gene_list.csv"

# Agent 3
copy_if_exists "${AGENT3_OUT}/ranked_pathway_list.csv"       "${FINAL}/ranked_pathway_list.csv"
copy_if_exists "${AGENT3_OUT}/cad_network_map.html"          "${FINAL}/cad_network_map.html"

# Agent 4
copy_if_exists "${AGENT4_OUT}/patient_subtype_assignments.csv" "${FINAL}/cad_subtype_model.csv"
copy_if_exists "${AGENT4_OUT}/subtype_outcome_associations.csv" "${FINAL}/subtype_outcome_associations.csv"

echo ""
echo "[consolidate] Final outputs:"
ls -lh "${FINAL}/"
echo ""
echo "Pull to Mac: bash scripts/sync_results.sh final"
