#!/bin/bash
#SBATCH -J rev1
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/rev1_%j.out
#SBATCH -e /path/to/cad-genetics/rev1_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
O=$_REPO/exports/agent1/mvmr_statinfix
PLINK=/path/to/software/plink
LDREF=$SCRATCH/lpa_ldref
mkdir -p $O $_REPO/exports/agent1/revision

echo "############ [1] MVMR: PAD (new, reviewer request) + CAD re-run for the instrument dump"
for T in PAD CAD; do
  echo "=== $T ==="
  Rscript $_REPO/scripts/agent1_genetics/R/run_mvmr.R --outcome $T --sumstat_dir $S \
    --out $O/mvmr_$T.csv
done

echo "############ [2] LD among the primary MVMR instruments"
# "independent by construction" is a distance rule, not an LD guarantee. Measure it.
IV=$O/mvmr_CAD.instruments.tsv
awk 'NR>1{print $1}' $IV > /tmp/mvmr_ivs.txt
wc -l < /tmp/mvmr_ivs.txt
$PLINK --bfile $LDREF --extract /tmp/mvmr_ivs.txt --r2 inter-chr --ld-window-r2 0 \
       --out /tmp/mvmr_ld --threads 16 > /dev/null 2>&1 || true
python3 - <<'PY'
import csv, os
ids = [l.split()[0] for l in open('/tmp/mvmr_ivs.txt')]
bim = set()
for l in open(os.environ.get('LDREF_BIM', '/path/to/scratch/lpa_ldref.bim')):
    bim.add(l.split('\t')[1])
present = [i for i in ids if i in bim]
print(f"[LD] instruments: {len(ids)}; present in LD reference: {len(present)}")
f = '/tmp/mvmr_ld.ld'
if os.path.exists(f):
    r2 = []
    for i, line in enumerate(open(f)):
        if i == 0: continue
        p = line.split()
        if p: r2.append(float(p[-1]))
    r2.sort()
    if r2:
        print(f"[LD] pairs={len(r2)} max_r2={r2[-1]:.5f} "
              f"p99={r2[int(0.99*len(r2))-1]:.5f} median={r2[len(r2)//2]:.6f}")
        print(f"[LD] pairs with r2>0.01: {sum(1 for x in r2 if x>0.01)}; r2>0.001: {sum(1 for x in r2 if x>0.001)}")
else:
    print("[LD] no .ld file produced")
PY

echo "############ [3] stroke-subtype contrasts + Lp(a) instrument sensitivity"
Rscript $_REPO/scripts/agent1_genetics/R/stroke_subtype_contrasts.R \
  --sumstat_dir $S --megastroke_dir $_REPO/data/external/megastroke \
  --ld_ref $LDREF --plink $PLINK --nboot 1000 \
  --out_prefix $_REPO/exports/agent1/revision/stroke_subtype

echo ALL_DONE
