#!/bin/bash
: "${DECODE_TOKEN:?}"
cd "$SCRATCH/pqtl2/decode"
while IFS=: read -r G FN; do
  [ -z "$G" ] && continue
  [ -s "$G.full.txt.gz" ] && { echo "$G present"; continue; }
  curl -sS --config - <<CFG
url = "https://download.decode.is/s3/download?token=${DECODE_TOKEN}&file=${FN}"
output = "${G}.full.txt.gz"
retry = 20
retry-delay = 5
max-time = 5400
CFG
  echo "$G $(stat -c%s "$G.full.txt.gz" 2>/dev/null || echo 0)" >> decode_dl.log
done < files.txt
echo DECODE_ALL_DONE >> decode_dl.log
