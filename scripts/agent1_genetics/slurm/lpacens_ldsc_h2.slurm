#!/bin/bash
#SBATCH -J lpacens_ldsc
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 01:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/lpacens_ldsc_%j.out
#SBATCH -e /path/to/cad-genetics/logs/lpacens_ldsc_%j.err
# Is the genome-wide inflation of the censored Lp(a) GWAS (lambdaGC 1.124 off the LPA locus, vs 1.007
# for the measured-only rerun) polygenic signal or confounding? LDSC intercept and ratio answer that.
# Three GWAS, same treatment: published (merged/), measured-only rerun (merged_lpameasctrl/), censored
# (merged_lpacens/). chr6:155-166 Mb is removed before munging, as in compare_genome_summary.txt.
# Univariate LDSC with the HapMap3 weight LD scores as the reference, the standard h2 set-up.
set -euo pipefail
umask 077
_REPO=/path/to/cad-genetics
source "${_REPO}/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_LDSC}"; set -u

GW="${AGENT1_OUT}/gwas_summary_stats"
OUT="${_REPO}/exports/agent1/lpa_censored/ldsc_h2"
WK="${SCRATCH}/lpacens_ldsc"
mkdir -p "$OUT" "$WK"
# The project's hm3 weights directory ships ldscore files only, with no .M/.M_5_50, so it cannot
# serve as --ref-ld-chr (job 3549578 died on exactly that). data/reference/ldsc_fixed/eur.* is a
# complete univariate EUR set, which is what an intercept and a total h2 need.
W="${_REPO}/data/reference/ldsc_fixed/eur."

for arm in published:merged measonly:merged_lpameasctrl censored:merged_lpacens; do
  name=${arm%%:*}; dir=${arm#*:}
  f="$GW/$dir/lipids_LPA.txt.gz"
  [ -s "$f" ] || { echo "MISSING $f"; exit 1; }
  echo "[$(date)] $name <- $f"
  # REGENIE columns: CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P EXTRA.
  # P from LOG10P, floored at 1e-300 so that munge does not drop the strongest rows as out of bounds.
  zcat "$f" | awk 'BEGIN{OFS="\t"; print "SNP","A1","A2","FRQ","N","Z","P"}
    NR>1 && $10!="NA" && $10>0 && !($1==6 && $2>=155000000 && $2<=166000000) {
      p=exp(-$12*log(10)); if (p<1e-300) p=1e-300
      print $3,$5,$4,$6,$7,$9/$10,p }' > "$WK/$name.txt"
  python "${MUNGE_PY}" --sumstats "$WK/$name.txt" --merge-alleles "${_REPO}/data/reference/ldsc_fixed/w_hm3.snplist" \
    --maf 0.01 --chunksize 500000 --out "$WK/$name"
  python "${LDSC_PY}" --h2 "$WK/$name.sumstats.gz" --ref-ld-chr "$W" --w-ld-chr "$W" --out "$OUT/h2_$name"
done

echo "== summary (chr6:155-166 Mb removed) =="
for name in published measonly censored; do
  printf "%-10s " "$name"
  grep -E "^Total Observed scale h2|^Lambda GC|^Mean Chi|^Intercept|^Ratio" "$OUT/h2_$name.log" | tr '\n' ' '; echo
done | tee "$OUT/summary.txt"
rm -f "$WK"/*.txt
echo ALL_DONE
