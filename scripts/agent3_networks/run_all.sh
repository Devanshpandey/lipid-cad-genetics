#!/usr/bin/env bash
# Agent 3 — Network Biology
# Prerequisite: Agent 2 complete (ranked_causal_gene_table.csv)
# TODO: implement network construction scripts

set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${REPO_DIR}/config/tacc_paths.sh"

REQUIRED="${AGENT2_OUT}/ranked_causal_gene_table.csv"
if [[ ! -f "${REQUIRED}" ]]; then
  echo "ERROR: Agent 2 output missing: ${REQUIRED}"
  echo "Run agent2_genes/run_all.sh first."
  exit 1
fi

echo "Agent 3 — Network Biology: coming soon"
echo "Will implement: STRING/BioGRID PPI, WGCNA co-expression,"
echo "  GSEA pathway enrichment, plaque scRNA-seq mapping"
