#!/usr/bin/env bash
# Push local code changes to TACC.
# Run from your Mac after editing scripts.
# Never pushes data/ or results/ — those stay on TACC.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/tacc_paths.sh"

echo "==> Syncing code to ${TACC_USER}@${TACC_HOST}:${TACC_WORK}/"

rsync -avz --progress \
  --exclude='.git/' \
  --exclude='results/' \
  --exclude='data/' \
  --exclude='logs/' \
  --exclude='*.DS_Store' \
  --exclude='__pycache__/' \
  --exclude='*.pyc' \
  --exclude='.Rhistory' \
  --exclude='.RData' \
  "${SCRIPT_DIR}/../" \
  "${TACC_USER}@${TACC_HOST}:${TACC_WORK}/"

echo "==> Done. Code is live on TACC."
echo "    Next: ssh ${TACC_USER}@${TACC_HOST}"
echo "    cd ${TACC_WORK}"
