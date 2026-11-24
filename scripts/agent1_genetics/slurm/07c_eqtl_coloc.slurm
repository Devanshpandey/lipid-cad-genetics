#!/bin/bash
#------------------------------------------------------------
# cis-eQTL colocalization for the lipid drug-target genes.
# Extends the earlier run (10 genes x 3 tissues x LDL-C only) to the full
# drug-target set, all five GTEx tissues, and the CAD outcome.
# Submit: sbatch 07c_eqtl_coloc.slurm
#------------------------------------------------------------
#SBATCH -J eqtl_coloc2
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 06:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -euo pipefail
_REPO="${SLURM_SUBMIT_DIR:-/path/to/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

REF="${_REPO}/data/reference/eqtl_catalogue"
WORK="${SCRATCH}/qtl_coloc"
mkdir -p "${WORK}"

# symbol  ENSG  chr(GRCh37)
cat > "${WORK}/genes.tsv" <<'G'
PCSK9	ENSG00000169174	1
SORT1	ENSG00000134243	1
CELSR2	ENSG00000143126	1
PSRC1	ENSG00000134222	1
ANGPTL3	ENSG00000132855	1
APOB	ENSG00000084674	2
HMGCR	ENSG00000113161	5
LPA	ENSG00000198670	6
NPC1L1	ENSG00000015520	7
TRIB1	ENSG00000173334	8
LPL	ENSG00000175445	8
APOC3	ENSG00000110245	11
APOA5	ENSG00000110243	11
LIPC	ENSG00000166035	15
CETP	ENSG00000087237	16
LDLR	ENSG00000130164	19
APOE	ENSG00000130203	19
APOC1	ENSG00000130208	19
G

TISSUES="Liver Artery_Aorta Artery_Coronary Heart_Left_Ventricle Whole_Blood"
ENS=$(awk '{print $2}' "${WORK}/genes.tsv" | paste -sd'|')

for TIS in ${TISSUES}; do
  SRC="${REF}/${TIS}.tsv.gz"
  DST="${WORK}/${TIS}_qtl.tsv.gz"
  if [ ! -f "${SRC}" ]; then echo "  MISSING ${SRC} — skipping"; continue; fi
  if [ -s "${DST}" ]; then echo "  ${TIS}: extract exists, reusing"; continue; fi
  echo "[$(date '+%H:%M:%S')] extracting ${TIS}"
  # Column 7 is gene_id in the eQTL Catalogue GTEx_V8/ge layout. Keep the
  # header so downstream naming is header-driven, never positional.
  zcat "${SRC}" | awk -F'\t' -v ids="${ENS}" \
    'BEGIN{n=split(ids,a,"|"); for(i=1;i<=n;i++) want[a[i]]=1}
     NR==1 || ($7 in want)' | gzip > "${DST}"
  echo "  ${TIS}: $(zcat "${DST}" | wc -l) rows, $(zcat "${DST}" | tail -n +2 | cut -f7 | sort -u | wc -l) genes found"
done

set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_R}"
set -u

Rscript "${SCRIPTS_DIR}/agent1_genetics/R/run_qtl_coloc.R" \
  --genes       "${WORK}/genes.tsv" \
  --qtl_dir     "${WORK}" \
  --tissues     "$(echo ${TISSUES} | tr ' ' ',')" \
  --sumstat_dir "${AGENT1_OUT}/gwas_summary_stats/merged" \
  --traits      "lipids_LDL_C,lipids_TRIGLY,lipids_HDL_C,lipids_APOB,outcomes_CAD" \
  --window_kb   500 \
  --out         "${WORK}/cis_eqtl_coloc_v2.csv"

echo "[$(date)] Done."
