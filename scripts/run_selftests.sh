#!/usr/bin/env bash
# Synthetic self-tests: no UK Biobank data are read. Run on a compute node.
#   bash scripts/run_selftests.sh [WORKDIR] [THREADS]
# Covers the genetic comparison with Lp(a)HORIZON and its post hoc analyses, the forecasts, the amplification
# test and three review analyses. Needs R with data.table, survival, optparse, splines, parallel, yaml, digest,
# MendelianRandomization and MVMR. Exit status 1 if any test fails.
set -uo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${1:-$(mktemp -d)}"; THREADS="${2:-4}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
A1="$HERE/scripts/agent1_genetics"; LH="$A1/lpa_horizon"; LS="$A1/lpa_susceptibility"; NR="$A1/nature_round"
mkdir -p "$WORK"; fails=0
run() {
  local name="$1"; shift
  local log="$WORK/$name.log" t0=$SECONDS
  if "$@" > "$log" 2>&1; then echo "PASS  $name  ($((SECONDS - t0)) s)"
  else echo "FAIL  $name  ($((SECONDS - t0)) s)  see $log"; tail -5 "$log" | sed 's/^/      /'; fails=$((fails + 1)); fi
}
for s in lpa_horizon lpa_horizon_posthoc lpa_horizon_posthoc2 lpa_forecast lpa_horizon_posthoc3; do
  run "$s" Rscript "$LH/$s.R" --step selftest --workdir "$(mktemp -d "$WORK/$s.XXXXXX")" --threads "$THREADS"
done
for s in lpa_susceptibility lpa_susceptibility_posthoc; do
  run "$s" Rscript "$LS/$s.R" --step selftest --workdir "$(mktemp -d "$WORK/$s.XXXXXX")" --threads "$THREADS"
done
run nat_apob_lpa_mvmr Rscript "$NR/nat_apob_lpa_mvmr_selftest.R" --workdir "$WORK/nat_apob_lpa" --threads "$THREADS"
run nat_statin_selection Rscript "$NR/nat_statin_selection.R" --selftest --selftest_dir "$WORK/nat_statin" --cores "$THREADS"
run nat_ldl_apob_robust Rscript "$NR/nat_ldl_apob_robust.R" --selftest --selftest_dir "$WORK/nat_robust" --ncores "$THREADS" --selftest_reps 10000
echo "$fails failed"
[ "$fails" -eq 0 ]
