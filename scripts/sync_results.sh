#!/usr/bin/env bash
# Pull results from TACC to your local Mac for visualization.
# Usage:
#   ./sync_results.sh              # pull all agents
#   ./sync_results.sh agent1       # pull only agent1_genetics
#   ./sync_results.sh agent1 agent2

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/tacc_paths.sh"

LOCAL_RESULTS="${SCRIPT_DIR}/../results"
mkdir -p "${LOCAL_RESULTS}"

AGENTS=("$@")
if [[ ${#AGENTS[@]} -eq 0 ]]; then
  AGENTS=("agent1_genetics" "agent2_genes" "agent3_networks" "agent4_subtypes" "final")
fi

for agent in "${AGENTS[@]}"; do
  # Normalize: accept "agent1" or "agent1_genetics"
  case "${agent}" in
    agent1|agent1_genetics) remote="agent1_genetics" ;;
    agent2|agent2_genes)    remote="agent2_genes" ;;
    agent3|agent3_networks) remote="agent3_networks" ;;
    agent4|agent4_subtypes) remote="agent4_subtypes" ;;
    final)                  remote="final" ;;
    *) echo "Unknown agent: ${agent}"; continue ;;
  esac

  echo "==> Pulling results/${remote} from TACC..."
  rsync -avz --progress \
    --exclude='*.loco' \
    --exclude='tmp_step1*' \
    "${TACC_USER}@${TACC_HOST}:${TACC_WORK}/results/${remote}/" \
    "${LOCAL_RESULTS}/${remote}/"
done

echo "==> Results synced to ${LOCAL_RESULTS}/"
