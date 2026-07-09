#!/usr/bin/env bash
# Agent 2 — Gene Prioritization
# Prerequisite: Agent 1 complete (finemapping_credible_sets.csv, coloc_results.csv)
#
# Usage:
#   nohup bash scripts/agent2_genes/run_all.sh \
#     > logs/agent2_$(date +%Y%m%d_%H%M%S).log 2>&1 &

set -euo pipefail
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
source "${REPO_DIR}/config/tacc_paths.sh"

REQUIRED="${AGENT1_OUT}/finemapping_credible_sets.csv"
if [[ ! -f "${REQUIRED}" ]]; then
  echo "ERROR: Agent 1 output missing: ${REQUIRED}"
  echo "Run agent1_genetics/run_all.sh first."
  exit 1
fi

echo "Agent 2 — Gene Prioritization starting..."
bash "${REPO_DIR}/scripts/agent2_genes/slurm/run_agent2.sh"
