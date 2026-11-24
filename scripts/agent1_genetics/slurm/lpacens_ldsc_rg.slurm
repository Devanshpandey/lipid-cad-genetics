#!/bin/bash
#SBATCH -J lpacens_rg
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 00:30:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/lpacens_rg_%j.out
#SBATCH -e /path/to/cad-genetics/logs/lpacens_rg_%j.err
# Follow-up to lpacens_ldsc_h2.slurm (job 3552601), which found the censored Lp(a) GWAS carries
# roughly three times the off-locus heritability of the measured-only rerun (0.0237 vs 0.0078) with
# an intercept of 1.0372. Heritability cannot rise from sample size alone, so the extra signal is
# either real trans-acting Lp(a) biology that the measured-only subsample loses to selection, or
# structure manufactured by the censored rank-normal coding (about 17% of people tie at the two
# assay boundaries, and who is out of range tracks LPA genotype).
#
# Genetic correlation separates those. Real biology correlates with the measured-only trait; coding
# artefact does not. Both arms share nearly all their samples (measured-only is a subset), so the rg
# intercept is expected to be large and is not interpretable here; rg itself is what matters, and
# LDSC estimates it in the presence of overlap. h2 for the measured-only arm is small, so expect a
# wide interval: a point estimate near 1 with a wide CI still argues for real signal, while a point
# estimate near 0 argues for artefact.
#
# Reuses the munged sumstats the h2 job left in $SCRATCH, so nothing is recomputed. chr6:155-166 Mb
# was already removed before those were munged.
set -euo pipefail
umask 077
_REPO=/path/to/cad-genetics
source "${_REPO}/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_LDSC}"; set -u

WK="${SCRATCH}/lpacens_ldsc"
OUT="${_REPO}/exports/agent1/lpa_censored/ldsc_h2"
W="${_REPO}/data/reference/ldsc_fixed/eur."
mkdir -p "$OUT"

for f in censored measonly published; do
  [ -s "$WK/$f.sumstats.gz" ] || { echo "MISSING $WK/$f.sumstats.gz: rerun lpacens_ldsc_h2.slurm first"; exit 1; }
done

python "${LDSC_PY}" --rg "$WK/censored.sumstats.gz,$WK/measonly.sumstats.gz,$WK/published.sumstats.gz" \
  --ref-ld-chr "$W" --w-ld-chr "$W" --out "$OUT/rg_censored_vs_measured"

echo "== rg summary (chr6:155-166 Mb removed; near-complete sample overlap, so the intercept is not a confounding read) =="
sed -n '/^Summary of Genetic Correlation Results/,/^$/p' "$OUT/rg_censored_vs_measured.log" | tee "$OUT/rg_summary.txt"
echo ALL_DONE
