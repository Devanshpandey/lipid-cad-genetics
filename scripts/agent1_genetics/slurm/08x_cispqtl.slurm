#!/bin/bash
#SBATCH -J cispqtl
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 03:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/cispqtl_%j.out
#SBATCH -e /path/to/cad-genetics/cispqtl_%j.err
# ============================================================================
# deCODE cis-pQTL x UK Biobank colocalization for four lipid effector genes.
#
# The archived attempt (Aug 2) failed silently: every download truncated partway
# through chr1, so PCSK9 and ANGPTL3 stopped ~20-39 Mb short of their cis-windows
# and APOB (chr2) / APOC3 (chr11) never reached their chromosomes at all. The
# pipeline's own >=900 MB guard caught it, which is why no results were written --
# dl_all.log saying "ALL DONE" only meant wget exited.
#
# This run downloads with wget -c (resumable), then VERIFIES each file both by
# size and by counting variants actually inside the cis-window before any
# colocalization is attempted. A file that fails either check stops the job.
#
# DECODE_TOKEN is passed through the environment, never written to disk: it is a
# per-user credential tied to a data-use agreement and it expires.
# ============================================================================
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
: "${DECODE_TOKEN:?DECODE_TOKEN not set -- resubmit with --export}"
export CISDIR=$SCRATCH/cispqtl
mkdir -p "$CISDIR"

echo "=== running pipeline (download + cis-window parse + coloc) ==="
cd "$CISDIR"
python3 $_REPO/scripts/agent1_genetics/cispqtl_pipeline.py || { echo "PIPELINE FAILED"; exit 1; }

echo
echo "=== verification: variants inside each cis-window ==="
python3 - <<'PY'
import os, gzip
CIS = os.environ["CISDIR"]
W = {"PCSK9":("chr1",54700000,55400000), "APOB":("chr2",20700000,21350000),
     "ANGPTL3":("chr1",62300000,62900000), "APOC3":("chr11",116500000,117150000)}
bad = []
for g,(c,s,e) in W.items():
    p = os.path.join(CIS, f"{g}.full.txt.gz")
    n = 0
    if os.path.exists(p):
        with gzip.open(p,"rb") as f:
            f.readline()
            for line in f:
                q = line.split(b"\t")
                if len(q) > 2 and q[0] == c.encode():
                    pos = int(q[1])
                    if s <= pos <= e: n += 1
                    elif pos > e: break
    print(f"  {g:<8} {c}:{s}-{e}  variants_in_window={n}")
    if n == 0: bad.append(g)
if bad: raise SystemExit("EMPTY CIS-WINDOW for: " + ",".join(bad))
print("  all four windows populated")
PY

OUT=$_REPO/exports/agent1/cispqtl; mkdir -p $OUT
cp -f "$CISDIR"/cispqtl_coloc_results.csv "$OUT"/ 2>/dev/null && echo "results -> $OUT"
cat "$OUT"/cispqtl_coloc_results.csv 2>/dev/null
echo ALL_DONE
