#!/bin/bash
#SBATCH -J pppcoloc
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/pppcoloc_%j.out
#SBATCH -e /path/to/cad-genetics/pppcoloc_%j.err
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
export PPPDIR=$SCRATCH/ukbppp
export OUTDIR=$SCRATCH/ukbppp
python3 $_REPO/scripts/agent1_genetics/ukbppp_coloc.py || { echo "FAILED"; exit 1; }
O=$_REPO/exports/agent1/cispqtl; mkdir -p $O
cp -f $SCRATCH/ukbppp/ukbppp_coloc_results.csv $O/ && echo "results -> $O"
echo; echo "=== SIDE BY SIDE: SomaScan (deCODE) vs Olink (UKB-PPP) ==="
python3 - <<'PY'
import pandas as pd, os
O=os.environ["_REPO"]+"/exports/agent1/cispqtl" if "_REPO" in os.environ else "/path/to/cad-genetics/exports/agent1/cispqtl"
d=pd.read_csv(O+"/cispqtl_coloc_results.csv"); p=pd.read_csv(O+"/ukbppp_coloc_results.csv")
d["cmp"]=d.comparison.str.replace("pQTL_vs_","",regex=False)
p["cmp"]=p.comparison.str.replace("pQTL_vs_","",regex=False)
m=d.merge(p,on=["gene","cmp"],suffixes=("_soma","_olink"),how="outer")
print(m[["gene","cmp","PP_H4_soma","PP_H4_olink","nsnps_soma","nsnps_olink"]].to_string(index=False))
PY
echo ALL_DONE
