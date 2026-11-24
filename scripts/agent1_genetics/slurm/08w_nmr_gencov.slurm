#!/bin/bash
#SBATCH -J nmrgc
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 01:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/nmrgc_%j.out
#SBATCH -e /path/to/cad-genetics/nmrgc_%j.err
# ============================================================================
# Table S28 (NMR LDL-P vs LDL-C) reports conditional F 3.4 / 3.3 and reads them as
# weak instruments. nmr_mvmr.R:48 computed them with gencov = 0 -- but LDL_P and
# LDL_C_NMR are derived from the SAME NMR spectrum in the SAME participants, so
# their SNP-exposure errors covary about as strongly as two measurements can.
# Correcting the identical defect elsewhere moved conditional F 7.5 -> 48 (primary
# ApoB model) and 4.8 -> 104 (statin-naive arm). If it does the same here, "weak
# instruments" is the wrong reading of Table S28 and the caption must change.
#
# The NMR summary statistics were absent from /work and recovered from the corral
# archive; the driver script was recovered from the same place.
# ============================================================================
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
ARCH=/path/to/archive/results/agent1_genetics
W=$SCRATCH/nmr_gencov; mkdir -p $W
OUT=$_REPO/exports/agent1/nmr_gencov; mkdir -p $OUT

echo "=== staging merged NMR sumstats from corral ==="
for T in LDL_P LDL_C_NMR; do
  SRC=$ARCH/gwas_summary_stats/nmr/nmr_${T}.regenie.gz
  [ -s "$SRC" ] || { echo "MISSING $SRC"; exit 1; }
  cp -n "$SRC" $W/ && echo "  staged nmr_${T}.regenie.gz ($(du -h $W/nmr_${T}.regenie.gz | cut -f1))"
  gunzip -c $W/nmr_${T}.regenie.gz > $W/${T}.txt
done

CAD=$_REPO/data/external/cardiogram/cad.add.160614.website.txt
PH=$_REPO/results/agent1_genetics/phenotypes/pheno_nmr.txt
[ -s "$CAD" ] || { echo "MISSING CARDIoGRAM"; exit 1; }
[ -s "$PH" ]  || { echo "MISSING pheno_nmr"; exit 1; }

echo "=== running NMR MVMR with and without the covariance correction ==="
Rscript $_REPO/scripts/agent1_genetics/R/nmr_mvmr.R \
  $W/LDL_P.txt $W/LDL_C_NMR.txt "$CAD" "$OUT" "$PH"
echo ALL_DONE
