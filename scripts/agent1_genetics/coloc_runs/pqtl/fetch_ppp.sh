#!/bin/bash
: "${SYNAPSE_TOKEN:?}"
cd "$SCRATCH/pqtl2/ukbppp"
while read -r G ID; do
  [ -z "$G" ] && continue
  [ -s "$G.tar" ] && { echo "$G present" >> ppp_dl.log; continue; }
  # header from stdin via --config so the bearer token stays out of argv
  curl -sSL --config - <<CFG
url = "https://repo-prod.prod.sagebase.org/repo/v1/entity/${ID}/file"
header = "Authorization: Bearer ${SYNAPSE_TOKEN}"
output = "${G}.tar"
retry = 10
max-time = 3600
CFG
  echo "$G $(stat -c%s "$G.tar" 2>/dev/null || echo 0)" >> ppp_dl.log
  mkdir -p "$G" && tar xf "$G.tar" -C "$G" 2>/dev/null && echo "  $G unpacked: $(find "$G" -type f | wc -l) files" >> ppp_dl.log
done < ents.txt
echo PPP_ALL_DONE >> ppp_dl.log
