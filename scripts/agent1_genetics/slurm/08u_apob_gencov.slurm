#!/bin/bash
#SBATCH -J apobgc
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 01:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/apobgc_%j.out
#SBATCH -e /path/to/cad-genetics/apobgc_%j.err
# ============================================================================
# Recompute conditional F and Q_A for the LDL-C/ApoB/triglyceride model with the
# SNP-exposure covariance supplied, instead of gencov = 0.
#
# WHY. All three exposure GWAS are the same 402,200 participants, so the
# SNP-exposure effect estimates covary. run_mvmr_apob.R:281 passed gencov = 0.
# mvmr_covariance.R exists to fix exactly this and its own header calls gencov=0
# a bug -- but its trait vector is c(LDL_C, HDL_C, TRIGLY, LPA) and ApoB was
# never added, so the ONE model whose headline claim is "conditional F < 10"
# is the one model that never got the correction. In the four-exposure model the
# correction raised conditional F by 32-100%; LDL-C and ApoB are far more
# correlated than any pair there, so the effect should be larger.
#
# WHAT THIS DOES NOT CHANGE. gencov enters strength_mvmr() and pleiotropy_mvmr()
# only. The mr_mvivw point estimates and standard errors -- and therefore the
# inversion (LDL-C OR 1.81 vs ApoB 0.72) that the non-identifiability claim rests
# on -- are untouched. This is a test of the REASON given, not of the finding.
#
# READ THE OUTPUT AS FOLLOWS.
#   cond_F_phenocov still < 10  -> the manuscript's "weak instruments" wording stands.
#   cond_F_phenocov >> 10       -> the wording is wrong and must become the stronger
#                                  claim: strongly instrumented, still unresolvable.
#   mvivw_b / mvivw_se moved    -> something else is wrong; stop and investigate.
# ============================================================================
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
PH=$_REPO/results/agent1_genetics/phenotypes/pheno_quantitative_statinfix.txt
OUT=$_REPO/exports/agent1/mvmr_apob_gencov
mkdir -p $OUT
ln -sf $_REPO/results/agent1_genetics/gwas_summary_stats/merged/lipids_APOB.txt.gz $S/lipids_APOB.txt.gz 2>/dev/null

for T in CAD MI MACE; do
  echo "########## $T ##########"
  Rscript $_REPO/scripts/agent1_genetics/R/run_mvmr_apob.R \
    --outcome $T --sumstat_dir $S --pheno $PH \
    --out $OUT/mvmr_apob_${T}.csv || echo "  $T FAILED"
done

echo
echo "########## COMPARISON vs the shipped uncorrected run ##########"
for T in CAD MI MACE; do
  OLD=$_REPO/exports/agent1/mvmr_apob_prunefix/mvmr_apob_${T}.csv
  NEW=$OUT/mvmr_apob_${T}.csv
  [ -f "$NEW" ] || continue
  python3 - "$OLD" "$NEW" "$T" <<'PY'
import csv, sys
old, new, T = sys.argv[1], sys.argv[2], sys.argv[3]
o = {r['exposure']: r for r in csv.DictReader(open(old))}
n = {r['exposure']: r for r in csv.DictReader(open(new))}
print(f"--- {T}")
for e in n:
    a, b = o.get(e), n[e]
    if not a: continue
    print(f"  {e:<7} cond_F {float(a['cond_F']):8.2f} -> phenocov {float(b['cond_F_phenocov']):8.2f}"
          f"   | mvivw_b {float(a['mvivw_b']):+.4f} -> {float(b['mvivw_b']):+.4f}"
          f"   se {float(a['mvivw_se']):.4f} -> {float(b['mvivw_se']):.4f}")
PY
done
echo ALL_DONE
