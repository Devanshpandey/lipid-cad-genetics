#!/bin/bash
#------------------------------------------------------------
# Fetch additional deCODE (SomaScan v4) cis-pQTL summary statistics.
#
# CREDENTIAL: DECODE_TOKEN must be supplied at submit time via
#   sbatch --export=ALL,DECODE_TOKEN=<token> 08aa_decode_fetch.slurm
# It is a per-user credential tied to the deCODE data-use agreement and it
# expires. Never write it into this file or into any output.
# The DUA also forbids redistributing these summary statistics, so keep the
# downloads on scratch and deposit only derived coloc statistics.
#------------------------------------------------------------
#SBATCH -J decode_fetch
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 10:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail
: "${DECODE_TOKEN:?DECODE_TOKEN not set — submit with --export=ALL,DECODE_TOKEN=...}"

OUT="${SCRATCH}/pqtl2/decode"
mkdir -p "${OUT}"
cd "${OUT}"

# gene:filename  (SomaScan seqid_sub_GENE_ProteinShortName)
FILES="
LDLR:13129_40_LDLR_LDLR.txt.gz
SORT1:11300_32_SORT1_SORT.txt.gz
APOE:2418_55_APOE_Apo_E.txt.gz
APOC1:15364_101_APOC1_Apo_C_I.txt.gz
APOA5:15363_32_APOA5_Apo_A_V.txt.gz
HMGCR:5230_99_HMGCR_HMGR.txt.gz
LDLRAP1:18183_3_LDLRAP1_ARH.txt.gz
APOA1:2750_3_APOA1_Apo_A_I.txt.gz
ANGPTL4:3796_79_ANGPTL4_ANGL4.txt.gz
LCAT:15413_3_LCAT_LCAT.txt.gz
"

echo "${FILES}" | grep -v '^$' | while IFS=: read -r GENE FN; do
  echo "${GENE} ${FN}"
done > manifest.txt

# -c is resumable; the URL carries the token so keep it out of logs (-nv would print it)
fetch_one() {
  GENE="$1"; FN="$2"
  DEST="${GENE}.full.txt.gz"
  if [ -s "${DEST}" ]; then echo "  ${GENE}: already present ($(stat -c%s "${DEST}") bytes)"; return 0; fi
  wget -c -q --tries=30 --timeout=90 --waitretry=5 \
       -O "${DEST}" "https://download.decode.is/s3/download?token=${DECODE_TOKEN}&file=${FN}"
  echo "  ${GENE}: $(stat -c%s "${DEST}" 2>/dev/null || echo 0) bytes"
}
export -f fetch_one
export DECODE_TOKEN

awk '{print $1, $2}' manifest.txt | xargs -P 4 -n 2 bash -c 'fetch_one "$0" "$1"'

echo "[$(date)] downloads complete:"
ls -la "${OUT}"/*.full.txt.gz | awk '{print "  ", $5, $9}'
