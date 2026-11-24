#!/bin/bash
#SBATCH -J natE_csusie
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/natE_coloc_susie_%j.out
#SBATCH -e /path/to/cad-genetics/logs/natE_coloc_susie_%j.err
# Nature round r1, Fig. 4a: coloc.susie at all five loci of the panel (published table: PCSK9 and LDLR only).
# Protocol, reproduction control and figure rule are in the header of nat_coloc_susie_5loci.R.
# Exit status 2 = the PCSK9/LDLR reproduction control failed and no new locus was run.
set -euo pipefail
umask 077
_REPO=/path/to/cad-genetics
source "$_REPO/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

SCRIPT=$_REPO/scripts/agent1_genetics/nature_round/nat_coloc_susie_5loci.R
OUT=$_REPO/exports/agent1/nature_round/coloc_susie5
for f in "$SCRIPT" "$_REPO/results/agent1_genetics/gwas_summary_stats/merged/lipids_LDL_C.txt.gz" \
         "$_REPO/results/agent1_genetics/gwas_summary_stats/merged/lipids_TRIGLY.txt.gz" \
         "$_REPO/results/agent1_genetics/gwas_summary_stats/merged/outcomes_CAD.txt.gz"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f"; exit 1; }
done
mkdir -p "$OUT"
echo "start $(date) on $(hostname)"
md5sum "$SCRIPT"
Rscript "$SCRIPT" "$OUT" "${SLURM_NTASKS:-16}"
echo "end $(date)"
