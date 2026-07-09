#!/usr/bin/env bash
# ============================================================
# Download all Agent 2 reference datasets to TACC
# All URLs verified 2026-05-23
#
# Usage:
#   nohup bash scripts/download_agent2_refs.sh \
#     > logs/download_agent2_$(date +%Y%m%d_%H%M%S).log 2>&1 &
# ============================================================

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${REPO_DIR}/config/tacc_paths.sh"

REF="${TACC_WORK}/data/reference"
mkdir -p "${REF}"/{gtex_v8,eqtl_catalogue,predictdb,pchic,opentargets,vep_cache,ukb_ppp}
mkdir -p "${LOGS_DIR}"

echo "=============================================="
echo " Agent 2 Reference Downloads"
echo " Destination: ${REF}"
echo " Started: $(date)"
echo "=============================================="

# ── 1. GTEx v8 eQTL — significant pairs (~3 GB) ───────────────
# New GCS bucket: adult-gtex (old gtex_analysis_v8 bucket is gone)
echo ""
echo "[$(date)] === 1. GTEx v8 significant eQTLs ==="
cd "${REF}/gtex_v8"
wget -c -q --show-progress \
  "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL.tar" \
  -O GTEx_Analysis_v8_eQTL.tar
echo "  Extracting..."
tar -xf GTEx_Analysis_v8_eQTL.tar --strip-components=1
echo "  Done: GTEx v8 significant eQTLs"

# Also grab the independent eQTLs (conditional analysis)
wget -c -q --show-progress \
  "https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/GTEx_Analysis_v8_eQTL_independent.tar" \
  -O GTEx_Analysis_v8_eQTL_independent.tar
tar -xf GTEx_Analysis_v8_eQTL_independent.tar --strip-components=1
echo "  Done: GTEx v8 independent eQTLs"

# ── 2. eQTL Catalogue — GTEx V8 all-pairs, tabix-indexed ─────
# GTEx all-pairs are NOT in GCS. Use eQTL Catalogue (harmonized,
# tabix-indexed by position — query specific loci without loading all).
# Key tissues for lipid/CAD biology:
#   Liver            — hepatic lipid metabolism (PCSK9, LDLR, APOB, SORT1)
#   Artery_Aorta     — vascular wall (PHACTR1, ADAMTS7)
#   Artery_Coronary  — coronary-specific (closest to CAD tissue)
#   Heart_Left_Ventricle — cardiac
#   Whole_Blood      — broad proxy, largest N in GTEx
echo ""
echo "[$(date)] === 2. eQTL Catalogue — GTEx V8 all-pairs (tabix, 5 tissues) ==="
BASE="https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/imported/GTEx_V8/ge"
for tissue in Liver Artery_Aorta Artery_Coronary Heart_Left_Ventricle Whole_Blood; do
  echo "  Downloading ${tissue} (~8-15 GB)..."
  wget -c -q --show-progress \
    "${BASE}/${tissue}.tsv.gz" \
    -O "${REF}/eqtl_catalogue/${tissue}.tsv.gz"
  wget -c -q --show-progress \
    "${BASE}/${tissue}.tsv.gz.tbi" \
    -O "${REF}/eqtl_catalogue/${tissue}.tsv.gz.tbi"
  echo "  Done: ${tissue}"
done
echo "  Done: eQTL Catalogue all-pairs"

# ── 3. PrediXcan MASHR-M eQTL models (~2 GB) ─────────────────
# For TWAS / S-PrediXcan. Note: URL uses /records/ not /record/
echo ""
echo "[$(date)] === 3. PrediXcan MASHR-M models ==="
cd "${REF}/predictdb"
wget -c -q --show-progress \
  "https://zenodo.org/records/3518299/files/mashr_eqtl.tar" \
  -O mashr_eqtl.tar
tar -xf mashr_eqtl.tar
echo "  Done: PrediXcan MASHR-M eQTL models"

# ── 4. Promoter-capture Hi-C — Javierre 2016 ─────────────────
# Blood cell type chromatin interactions from GEO GSE86731
# RAW.tar contains the peak matrix and interaction files
echo ""
echo "[$(date)] === 4. PCHi-C Javierre 2016 ==="
cd "${REF}/pchic"
wget -c -q --show-progress \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE86nnn/GSE86731/suppl/GSE86731_RAW.tar" \
  -O GSE86731_RAW.tar
tar -xf GSE86731_RAW.tar
# Extract just the peak matrix (the key file for interaction scoring)
gunzip -f PCHiC_peak_matrix_cutoff5.tsv.gz 2>/dev/null || \
  find . -name "*peak_matrix*" | head -1 | xargs -I{} gunzip -f {} || true
echo "  Done: PCHi-C"

# ── 5. OpenTargets Genetics — L2G scores (parquet) ───────────
# Locus-to-gene ML scores integrating eQTL, PCHi-C, distance, VEP
echo ""
echo "[$(date)] === 5. OpenTargets Genetics L2G ==="
mkdir -p "${REF}/opentargets/l2g"
# Download all parquet shards (multiple files, ~25 GB total)
OT_L2G="https://ftp.ebi.ac.uk/pub/databases/opentargets/genetics/latest/l2g"
# Get file list and download each parquet
wget -q -O - "${OT_L2G}/" 2>/dev/null | \
  grep -oP '(?<=href=")[^"]+\.parquet' | \
  while read f; do
    wget -c -q --show-progress "${OT_L2G}/${f}" -O "${REF}/opentargets/l2g/${f}"
  done || \
  # Fallback: direct wget recursive
  wget -r -np -nH --cut-dirs=8 -q --show-progress \
    --accept="*.parquet" \
    "${OT_L2G}/" -P "${REF}/opentargets/l2g/"
echo "  Done: OpenTargets L2G"

# ── 6. OpenTargets Platform 25.12 — targets + drugs ──────────
# Targets: tractability, gene info
# known_drug: approved drugs per target
# l2g_prediction: GWAS-to-gene predictions
echo ""
echo "[$(date)] === 6. OpenTargets Platform (targets, drugs) ==="
OT_PLAT="https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output"
for dtype in target known_drug l2g_prediction; do
  mkdir -p "${REF}/opentargets/${dtype}"
  echo "  Downloading ${dtype}..."
  wget -r -np -nH --cut-dirs=7 -q --show-progress \
    --accept="*.json.gz,*.parquet" \
    "${OT_PLAT}/${dtype}/" \
    -P "${REF}/opentargets/${dtype}/" || \
  rsync -rtz \
    "rsync://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/${dtype}/" \
    "${REF}/opentargets/${dtype}/" 2>/dev/null || true
  echo "  Done: ${dtype}"
done

# ── 7. Ensembl VEP cache GRCh38 v114 (~14 GB) ────────────────
echo ""
echo "[$(date)] === 7. Ensembl VEP cache GRCh38 ==="
cd "${REF}/vep_cache"
wget -c -q --show-progress \
  "https://ftp.ensembl.org/pub/release-114/variation/indexed_vep_cache/homo_sapiens_vep_114_GRCh38.tar.gz" \
  -O homo_sapiens_vep_114_GRCh38.tar.gz
echo "  Extracting VEP cache (~14 GB, takes a few minutes)..."
tar -xzf homo_sapiens_vep_114_GRCh38.tar.gz
echo "  Done: VEP cache"

# ── 8. UKB-PPP pQTL — sample + instructions for bulk ─────────
# Full dataset: ~80 GB, 1463 proteins across GCST90271001-GCST90274093
# Bulk download requires looping; we write a helper script here
echo ""
echo "[$(date)] === 8. UKB-PPP pQTL (writing bulk-fetch helper) ==="
mkdir -p "${REF}/ukb_ppp"
cat > "${REF}/ukb_ppp/fetch_ukb_ppp.sh" << 'PQTL_EOF'
#!/usr/bin/env bash
# Bulk fetch UKB-PPP pQTL summary stats (Sun et al. 2023, Nature)
# ~1463 proteins, GCST90271001-GCST90274093 range on EBI GWAS Catalog
# Run this separately: nohup bash fetch_ukb_ppp.sh > fetch_ukb_ppp.log 2>&1 &
set -euo pipefail
OUT_DIR="$(dirname "$0")"
BASE="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics"
for gcst_num in $(seq 90271001 90274093); do
  gcst="GCST${gcst_num}"
  range="GCST$(( (gcst_num / 1000) * 1000 + 1 ))-GCST$(( (gcst_num / 1000 + 1) * 1000 ))"
  url="${BASE}/${range}/${gcst}/harmonised/${gcst}.h.tsv.gz"
  out="${OUT_DIR}/${gcst}.h.tsv.gz"
  [[ -f "${out}" ]] && continue
  wget -q -c "${url}" -O "${out}" 2>/dev/null || true
done
PQTL_EOF
chmod +x "${REF}/ukb_ppp/fetch_ukb_ppp.sh"
echo "  UKB-PPP helper written to: ${REF}/ukb_ppp/fetch_ukb_ppp.sh"
echo "  Run separately (80 GB): nohup bash ${REF}/ukb_ppp/fetch_ukb_ppp.sh &"

# ── Summary ───────────────────────────────────────────────────
echo ""
echo "=============================================="
echo " Download complete: $(date)"
echo " Storage used:"
du -sh "${REF}"/*/  2>/dev/null
echo "=============================================="
