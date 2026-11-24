#!/bin/bash
#------------------------------------------------------------
# Build the LD panel used to define colocalization sentinels.
#
# PLINK --clump only ever considers variants that appear in the association
# file, so the panel needs exactly the union of genome-wide-significant
# variants across the lipid exposures -- not all 14.8M. Extracting them once
# turns a 1.4 TB random-read on corral into a ~15 GB sequential file on
# scratch, and (important) guarantees full panel coverage: --clump silently
# DROPS assoc variants missing from the panel, which is how earlier cis-MR
# instrument counts swung by panel density.
#
# Submit: sbatch 07a_coloc_ldpanel.slurm
#------------------------------------------------------------
#SBATCH -J coloc_ldpanel
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

MERGED="${AGENT1_OUT}/gwas_summary_stats/merged"
WORKDIR="${SCRATCH}/coloc_rerun"
mkdir -p "${WORKDIR}"

LIPIDS=(LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB)
# -log10(5e-8)
THRESH=7.30103

echo "[$(date '+%H:%M:%S')] Extracting GWS variant IDs per exposure"
: > "${WORKDIR}/gws_union.txt"
for TR in "${LIPIDS[@]}"; do
  zcat "${MERGED}/lipids_${TR}.txt.gz" \
    | awk -v t="${THRESH}" 'NR>1 && $12>t {print $3}' \
    > "${WORKDIR}/gws_${TR}.ids"
  echo "  ${TR}: $(wc -l < "${WORKDIR}/gws_${TR}.ids") GWS variants"
  cat "${WORKDIR}/gws_${TR}.ids" >> "${WORKDIR}/gws_union.txt"
done

sort -u "${WORKDIR}/gws_union.txt" -o "${WORKDIR}/gws_union.txt"
echo "[$(date '+%H:%M:%S')] Union: $(wc -l < "${WORKDIR}/gws_union.txt") unique variants"

echo "[$(date '+%H:%M:%S')] Extracting panel (all ${?} samples retained)"
"${PLINK2}" \
  --bfile "${UKB_GENO_STEP2}" \
  --extract "${WORKDIR}/gws_union.txt" \
  --make-bed \
  --out "${WORKDIR}/coloc_ldref" \
  --threads 8 \
  --memory 48000

echo "[$(date '+%H:%M:%S')] Panel built:"
wc -l < "${WORKDIR}/coloc_ldref.bim" | xargs echo "  variants:"
wc -l < "${WORKDIR}/coloc_ldref.fam" | xargs echo "  samples: "

# Coverage check: every GWS variant must be present, or --clump would silently
# undercount sentinels exactly the way the cis-MR panels did.
MISSING=$(comm -23 "${WORKDIR}/gws_union.txt" <(cut -f2 "${WORKDIR}/coloc_ldref.bim" | sort -u) | wc -l)
echo "  GWS variants absent from panel: ${MISSING}"
if [ "${MISSING}" -gt 0 ]; then
  echo "  NOTE: absent variants cannot index a clump; they are listed in gws_absent.txt"
  comm -23 "${WORKDIR}/gws_union.txt" <(cut -f2 "${WORKDIR}/coloc_ldref.bim" | sort -u) > "${WORKDIR}/gws_absent.txt"
fi

echo "[$(date)] Done."
