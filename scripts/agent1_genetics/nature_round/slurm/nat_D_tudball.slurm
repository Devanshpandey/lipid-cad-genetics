#!/bin/bash
#SBATCH -J nat_D_tudball
#SBATCH -p vm-small
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 08:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/nat_D_tudball_%j.out
#SBATCH -e /path/to/cad-genetics/logs/nat_D_tudball_%j.err
# Nature round, analysis D: Tudball et al. 2026 design on the paper's scores.
# How results are read: scripts/agent1_genetics/nature_round/PREREGISTERED_READOUT.md section D.
#
# The R script stops with status 3 if the port's published rows or the paper's Table S35 do not reproduce.
# That failure must end this job visibly, hence set -e and no "|| echo FAILED".
# Output is summary-level only; individual-level inputs stay on /scratch and /corral.
set -euo pipefail
umask 077
_REPO=/path/to/cad-genetics
source "$_REPO/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
# glmnet and sandwich live here, and this is the library the port ran against: PC1's IPTW rows need this build.
export R_LIBS_USER=/path/to/R/library
# The port ran with whatever BLAS threads a whole node gave it. glm uses LINPACK and glmnet its own Fortran,
# so this does not touch the estimates; it keeps a lone OpenBLAS pool from oversubscribing the VM slice.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

COHORT=/path/to/scratch/tudball-port/cohort_with_outcomes.csv.gz
W=/path/to/scratch/statint
RAW=/path/to/ukb_phenotypes
COVAR=$_REPO/results/agent1_genetics/phenotypes/covariates.txt
OUT=$_REPO/exports/agent1/nature_round/tudball
SCRIPT=$_REPO/scripts/agent1_genetics/nature_round/nat_tudball_paper_scores.R

for f in "$COHORT" "$W/dose_LPA.raw" "$W/dose_LDL.raw" "$W/lpa_score.txt" "$W/ldl_score.txt" "$COVAR" "$SCRIPT" \
         "$RAW/fid53.csv" "$RAW/fid40000.csv" "$RAW/fid6153.csv" "$RAW/fid6177.csv" \
         "$RAW/fid131296.csv" "$RAW/fid131298.csv" "$RAW/fid131300.csv" "$RAW/fid131302.csv" "$RAW/fid131304.csv" "$RAW/fid131306.csv"; do
  [ -s "$f" ] || { echo "MISSING INPUT: $f" >&2; exit 1; }
done
# Optional arm; the R script reports NA with a note if absent, so this only warns.
for f in "$RAW/fid30790.csv" "$RAW/fid30796.csv"; do
  [ -s "$f" ] || echo "NOTE: $f absent, repeat Lp(a) arm will report not_run" >&2
done
Rscript -e 'for (p in c("data.table","optparse","survival","glmnet","sandwich")) if (!requireNamespace(p, quietly = TRUE)) stop("missing R package: ", p); cat(R.version.string, " glmnet", as.character(packageVersion("glmnet")), " survival", as.character(packageVersion("survival")), " sandwich", as.character(packageVersion("sandwich")), "\n")'
mkdir -p "$OUT"
Rscript "$SCRIPT" \
  --cohort  "$COHORT" \
  --scores  "$W/dose_LPA.raw,$W/dose_LDL.raw" \
  --weights "$W/lpa_score.txt,$W/ldl_score.txt" \
  --raw_dir "$RAW" \
  --covar   "$COVAR" \
  --outdir  "$OUT"
echo ALL_DONE
