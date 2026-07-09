#!/usr/bin/env bash
# Agent 4 — Patient Subtyping
# Prerequisite: Agents 1-3 complete + UKB individual-level phenotype access
# TODO: implement subtyping scripts

set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${REPO_DIR}/config/tacc_paths.sh"

for f in \
  "${AGENT1_OUT}/prs_weights" \
  "${AGENT2_OUT}/ranked_causal_gene_table.csv" \
  "${AGENT3_OUT}/network_modules.json"; do
  if [[ ! -e "${f}" ]]; then
    echo "ERROR: Missing prerequisite: ${f}"
    exit 1
  fi
done

echo "Agent 4 — Patient Subtyping: coming soon"
echo "Will implement: pathway-specific PRS, LCA, NMF, Cox PH outcome associations"
