#!/bin/bash
#SBATCH -J nmr_sep
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 2:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
# Batch equivalent of run_nmr_separation.sh: merge -> LDSC rg -> two-trait MVMR.
set -uo pipefail
_REPO="/path/to/cad-genetics"
source "${_REPO}/config/tacc_paths.sh"
source "${CONDA_BASE}/etc/profile.d/conda.sh"
NMR="${AGENT1_OUT}/gwas_summary_stats/nmr"
OUT="${AGENT1_OUT}/strengthen/nmr_separation"; mkdir -p "$OUT"
WP="${LDSC_WEIGHTS_DIR}/weights.hm3_noMHC."
[ -f "${WP}1.l2.ldscore.gz" ] || WP="${LDSC_REF_PREFIX}"
MA=""; [ -f "${LDSC_HAPMAP3_SNPS}" ] && MA="--merge-alleles ${LDSC_HAPMAP3_SNPS}"

echo "===== 1) MERGE per-chr -> genome-wide ====="
for T in LDL_P LDL_C_NMR; do
  M="${NMR}/nmr_${T}.regenie.gz"
  if [ ! -s "$M" ]; then
    { zcat "${NMR}/nmr_chr1_${T}.regenie.gz" | head -1
      for c in $(seq 1 22); do zcat "${NMR}/nmr_chr${c}_${T}.regenie.gz" | tail -n +2; done; } | gzip > "$M"
  fi
  echo "[merge] ${T}: $(zcat "$M" | wc -l) rows"
done

echo "===== 2) MUNGE + LDSC genetic correlation ====="
conda activate ldsc
for T in LDL_P LDL_C_NMR; do
  zcat "${NMR}/nmr_${T}.regenie.gz" | awk 'NR==1{print "SNP\tA1\tA2\tN\tBETA\tP"; next}{printf "%s\t%s\t%s\t%s\t%s\t%g\n",$3,$5,$4,$7,$9,10^(-$12)}' > "${OUT}/munge_${T}.txt"
  python "${MUNGE_PY}" --sumstats "${OUT}/munge_${T}.txt" --snp SNP --a1 A1 --a2 A2 --N-col N \
    --signed-sumstats BETA,0 --p P ${MA} --out "${OUT}/${T}" 2>&1 | tail -3
done
python "${LDSC_PY}" --rg "${OUT}/LDL_P.sumstats.gz,${OUT}/LDL_C_NMR.sumstats.gz" \
  --ref-ld-chr "${LDSC_REF_PREFIX}" --w-ld-chr "${WP}" --out "${OUT}/nmr_rg" 2>&1 | tail -4
grep -A4 'Summary of Genetic Correlation' "${OUT}/nmr_rg.log" 2>/dev/null || true

echo "===== 3) MVMR: LDL_P vs LDL_C_NMR -> CARDIoGRAM CAD ====="
conda activate "${CONDA_ENV_R}"
Rscript "${_REPO}/scripts/agent1_genetics/strengthen/nmr_mvmr.R" \
  "${NMR}/nmr_LDL_P.regenie.gz" "${NMR}/nmr_LDL_C_NMR.regenie.gz" \
  "${_REPO}/data/external/cardiogram/cad.add.160614.website.txt" "${OUT}"
echo "===== DONE ====="; ls -la "${OUT}"
