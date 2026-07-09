#!/bin/bash
#------------------------------------------------------------
# prepare_burden_masks.sh
#
# One-time setup: build REGENIE annotation files for burden testing
# from UKB WES variant annotations.
#
# Source priority (checked in order):
#   1. joyce_annotated_vcf/  — VEP-annotated VCFs already on corral-repl (best)
#   2. wes_allchr.pvar       — parse from pvar directly (fallback, less annotation)
#
# Outputs (in ${AGENT1_OUT}/burden_masks/):
#   ukb_wes_anno.txt     — REGENIE --anno-file (CHROM POS REF ALT GENE ANNOTATION)
#   ukb_wes_setlist.txt  — REGENIE --set-list  (GENE CHROM START END VARIANT_LIST)
#   ukb_wes_masks.txt    — REGENIE --mask-def  (MASK_NAME ANNOTATION_CATS)
#
# Annotation categories written:
#   pLoF      — LOFTEE HC stop_gain / frameshift / essential splice
#   missense3 — damaging missense (REVEL ≥ 0.7 OR CADD ≥ 25, equivalent to MiS 0.7-1.0)
#   missense2 — moderate missense (REVEL 0.5-0.7 OR CADD 20-25)
#   synonymous — synonymous (negative control)
#
# Mask definitions:
#   mask1 = pLoF only
#   mask2 = pLoF + missense3  (as used in Koyama Nat Genet 2026)
#   mask3 = pLoF + missense3 + missense2
#
# Usage (run on login node or via srun — takes ~30-90 min for WES):
#   bash scripts/agent1_genetics/slurm/prepare_burden_masks.sh
#------------------------------------------------------------

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../config/tacc_paths.sh"

set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_PY}"
set -u

OUT_DIR="${UKB_BURDEN_ANNO_DIR}"   # set in tacc_paths.sh → ${AGENT1_OUT}/burden_masks
mkdir -p "${OUT_DIR}" "${LOGS_DIR}"

LOG_FILE="${LOGS_DIR}/prepare_burden_masks_$(date +%Y%m%d).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

ANNO_OUT="${UKB_BURDEN_ANNO_FILE}"
SETLIST_OUT="${UKB_BURDEN_SETLIST}"
MASKS_OUT="${UKB_BURDEN_MASKS}"

echo "[$(date)] Preparing REGENIE burden annotation files"
echo "  WES pfile:        ${UKB_WES_PFILE}.pgen/pvar/psam"
echo "  Joyce VCF dir:    ${UKB_WES_JOYCE_VCF}/"
echo "  Output dir:       ${OUT_DIR}/"

# ============================================================
# 1. Check joyce_annotated_vcf/ — preferred source
# ============================================================
JOYCE_DIR="${UKB_WES_JOYCE_VCF}"

if [[ -d "${JOYCE_DIR}" ]]; then
  echo "[$(date)] Found joyce_annotated_vcf/ — inspecting contents..."
  ls -lh "${JOYCE_DIR}" | head -30

  # Detect what files are present
  N_VCF=$(find "${JOYCE_DIR}" -name "*.vcf.gz" -o -name "*.vcf" | wc -l)
  N_TSV=$(find "${JOYCE_DIR}" -name "*.tsv*" -o -name "*.txt*" | wc -l)
  echo "  VCF/VCF.gz files: ${N_VCF}"
  echo "  TSV/TXT files:    ${N_TSV}"

  # Peek at first VCF or TSV header to understand format
  FIRST_VCF=$(find "${JOYCE_DIR}" -name "*.vcf.gz" | sort | head -1)
  FIRST_TSV=$(find "${JOYCE_DIR}" -name "*.tsv.gz" -o -name "*.tsv" | sort | head -1)

  USE_JOYCE_VCF=0
  USE_JOYCE_TSV=0

  if [[ -n "${FIRST_VCF}" ]]; then
    echo "  Peeking at: ${FIRST_VCF}"
    # Check VCF INFO fields for LOFTEE (LoF) and REVEL annotations
    if "${BCFTOOLS}" view -h "${FIRST_VCF}" 2>/dev/null | grep -qiE "LoF|LOFTEE|REVEL|CSQ|VEP"; then
      echo "  -> VCF has LOFTEE/REVEL annotations — using Joyce VCF as annotation source"
      USE_JOYCE_VCF=1
    else
      echo "  -> VCF INFO fields do not contain LOFTEE/REVEL — checking for TSV"
    fi
  fi

  if [[ -n "${FIRST_TSV}" && "${USE_JOYCE_VCF}" -eq 0 ]]; then
    echo "  Peeking at: ${FIRST_TSV}"
    USE_JOYCE_TSV=1
  fi

  if [[ "${USE_JOYCE_VCF}" -eq 1 ]]; then
    echo "[$(date)] Parsing Joyce annotated VCFs (all chromosomes)..."
    python - "${JOYCE_DIR}" "${ANNO_OUT}" "${SETLIST_OUT}" << 'PYEOF'
import sys, os, gzip, re
from collections import defaultdict

joyce_dir   = sys.argv[1]
anno_out    = sys.argv[2]
setlist_out = sys.argv[3]

# Find all VCF files, sort by chromosome
vcf_files = sorted(
    [os.path.join(joyce_dir, f) for f in os.listdir(joyce_dir)
     if f.endswith(".vcf.gz") or f.endswith(".vcf")],
    key=lambda x: (int(re.search(r"chr?(\d+)", os.path.basename(x)).group(1))
                   if re.search(r"chr?(\d+)", os.path.basename(x)) else 99)
)
if not vcf_files:
    vcf_files = sorted(
        [os.path.join(joyce_dir, f) for f in os.listdir(joyce_dir)
         if f.endswith(".gz") and "ann" in f.lower()],
    )
print(f"  Found {len(vcf_files)} VCF files")

# Parse CSQ/INFO field headers from first file to identify column positions
def get_csq_fields(vcf_path):
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("##INFO=<ID=CSQ") or line.startswith("##INFO=<ID=ANN"):
                m = re.search(r"Format: ([^\"]+)\"", line)
                if m:
                    return [f.strip() for f in m.group(1).split("|")]
            if not line.startswith("#"):
                break
    return []

def classify_variant(gene, consequence, lof, lof_filter, revel, cadd):
    # Priority: pLoF > missense3 > missense2 > synonymous
    if lof == "HC" and (not lof_filter or lof_filter == "."):
        return "pLoF"
    if any(c in consequence for c in
           ["stop_gained","frameshift","splice_donor_variant","splice_acceptor_variant",
            "start_lost","stop_lost"]):
        if lof in ("", ".", "LC"):
            # Check if splice AI-like or conserved — treat as pLoF if no filter
            if "splice" in consequence and not lof_filter:
                return "pLoF"
    if "missense" in consequence or "protein_altering" in consequence:
        r = float(revel) if revel and revel != "." else 0.0
        c = float(cadd)  if cadd  and cadd  != "." else 0.0
        if r >= 0.7 or c >= 25:  return "missense3"
        if r >= 0.5 or c >= 20:  return "missense2"
        return None  # low-confidence missense excluded
    if "synonymous" in consequence:
        return "synonymous"
    return None

anno_rows  = []
gene_vars  = defaultdict(list)  # gene → [(chrom, pos, ref, alt)]
gene_coords = defaultdict(lambda: {"chrom": "", "min": float("inf"), "max": 0})

for vcf_path in vcf_files:
    fname = os.path.basename(vcf_path)
    print(f"  Parsing {fname}...", flush=True)
    csq_fields = get_csq_fields(vcf_path)
    gene_idx = csq_fields.index("Gene") if "Gene" in csq_fields else \
               csq_fields.index("SYMBOL") if "SYMBOL" in csq_fields else 3
    cons_idx  = csq_fields.index("Consequence") if "Consequence" in csq_fields else 1
    lof_idx   = csq_fields.index("LoF") if "LoF" in csq_fields else -1
    loff_idx  = csq_fields.index("LoF_filter") if "LoF_filter" in csq_fields else -1
    rev_idx   = csq_fields.index("REVEL_score") if "REVEL_score" in csq_fields else -1
    cadd_idx  = csq_fields.index("CADD_phred") if "CADD_phred" in csq_fields else -1

    opener = gzip.open if vcf_path.endswith(".gz") else open
    n_written = 0
    with opener(vcf_path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            chrom = chrom.lstrip("chr")
            # Skip multi-allelics (already split by VEP typically)
            if "," in alt:
                alt = alt.split(",")[0]
            info = parts[7]
            csq_entry = ""
            for tok in info.split(";"):
                if tok.startswith("CSQ=") or tok.startswith("ANN="):
                    csq_entry = tok.split("=", 1)[1].split(",")[0]
                    break
            if not csq_entry or not csq_fields:
                continue
            csq_parts = csq_entry.split("|")
            def get_field(idx):
                return csq_parts[idx] if idx >= 0 and idx < len(csq_parts) else ""
            gene       = get_field(gene_idx)
            consequence = get_field(cons_idx)
            lof         = get_field(lof_idx)
            lof_filter  = get_field(loff_idx)
            revel       = get_field(rev_idx)
            cadd        = get_field(cadd_idx)
            if not gene or gene == ".":
                continue
            anno = classify_variant(gene, consequence, lof, lof_filter, revel, cadd)
            if anno is None:
                continue
            anno_rows.append(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{anno}")
            var_id = f"{chrom}:{pos}:{ref}:{alt}"
            gene_vars[gene].append(var_id)
            p = int(pos)
            gene_coords[gene]["chrom"] = chrom
            gene_coords[gene]["min"]   = min(gene_coords[gene]["min"], p)
            gene_coords[gene]["max"]   = max(gene_coords[gene]["max"], p)
            n_written += 1
    print(f"    {n_written:,} variants written from {fname}", flush=True)

print(f"\n  Total: {len(anno_rows):,} annotated variants, {len(gene_vars):,} genes")

# Write anno file
with open(anno_out, "w") as fh:
    fh.write("\n".join(anno_rows) + "\n")
print(f"  Anno file: {anno_out}")

# Write set list
with open(setlist_out, "w") as fh:
    for gene, vars_ in sorted(gene_vars.items()):
        coords = gene_coords[gene]
        var_str = ",".join(vars_)
        fh.write(f"{gene}\t{coords['chrom']}\t{coords['min']}\t{coords['max']}\t{var_str}\n")
print(f"  Set list: {setlist_out}")
PYEOF
    echo "[$(date)] Joyce VCF parsing complete."

  elif [[ "${USE_JOYCE_TSV}" -eq 1 ]]; then
    echo "[$(date)] Parsing Joyce TSV annotation file..."
    python - "${FIRST_TSV}" "${ANNO_OUT}" "${SETLIST_OUT}" << 'PYEOF'
import sys, gzip, os
from collections import defaultdict

tsv_path    = sys.argv[1]
anno_out    = sys.argv[2]
setlist_out = sys.argv[3]

opener = gzip.open if tsv_path.endswith(".gz") else open
print(f"  Reading: {tsv_path}")

with opener(tsv_path, "rt") as fh:
    header = fh.readline().rstrip().split("\t")
    print(f"  Columns: {header}")

# Detect column names (case-insensitive)
def find_col(header, candidates):
    for c in candidates:
        for i, h in enumerate(header):
            if h.lower() == c.lower():
                return i
    return -1

with opener(tsv_path, "rt") as fh:
    header = fh.readline().rstrip().split("\t")
    ci = {
        "chrom": find_col(header, ["CHROM","CHR","#CHROM"]),
        "pos":   find_col(header, ["POS","POSITION","BP"]),
        "ref":   find_col(header, ["REF","A2","ALLELE0"]),
        "alt":   find_col(header, ["ALT","A1","ALLELE1"]),
        "gene":  find_col(header, ["GENE","SYMBOL","GENENAME"]),
        "lof":   find_col(header, ["LOF","LoF","LOFTEE"]),
        "revel": find_col(header, ["REVEL","REVEL_score"]),
        "cadd":  find_col(header, ["CADD","CADD_phred"]),
        "cons":  find_col(header, ["consequence","CONSEQUENCE","Consequence"]),
        "anno":  find_col(header, ["annotation","ANNOTATION","category"]),
    }
    print(f"  Column map: {ci}")

    anno_rows  = []
    gene_vars  = defaultdict(list)
    gene_coords = defaultdict(lambda: {"chrom": "", "min": float("inf"), "max": 0})

    for line in fh:
        p = line.rstrip().split("\t")
        def g(k): return p[ci[k]].strip() if ci[k] >= 0 and ci[k] < len(p) else ""
        chrom = g("chrom").lstrip("chr")
        pos   = g("pos")
        ref   = g("ref")
        alt   = g("alt")
        gene  = g("gene")
        if not all([chrom, pos, ref, alt, gene]):
            continue
        # Use pre-computed annotation if available, otherwise derive
        anno = g("anno")
        if not anno:
            lof   = g("lof")
            cons  = g("cons")
            revel = g("revel")
            cadd  = g("cadd")
            if lof == "HC":
                anno = "pLoF"
            elif "missense" in cons:
                r = float(revel) if revel and revel != "." else 0.0
                c = float(cadd)  if cadd  and cadd  != "." else 0.0
                if r >= 0.7 or c >= 25: anno = "missense3"
                elif r >= 0.5 or c >= 20: anno = "missense2"
                else: continue
            elif "synonymous" in cons:
                anno = "synonymous"
            else:
                continue
        anno_rows.append(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{anno}")
        var_id = f"{chrom}:{pos}:{ref}:{alt}"
        gene_vars[gene].append(var_id)
        pp = int(pos)
        gene_coords[gene]["chrom"] = chrom
        gene_coords[gene]["min"]   = min(gene_coords[gene]["min"], pp)
        gene_coords[gene]["max"]   = max(gene_coords[gene]["max"], pp)

print(f"  {len(anno_rows):,} variants, {len(gene_vars):,} genes")
with open(anno_out, "w") as fh:
    fh.write("\n".join(anno_rows) + "\n")
with open(setlist_out, "w") as fh:
    for gene, vars_ in sorted(gene_vars.items()):
        coords = gene_coords[gene]
        fh.write(f"{gene}\t{coords['chrom']}\t{coords['min']}\t{coords['max']}\t{','.join(vars_)}\n")
print(f"  Written: {anno_out}")
print(f"  Written: {setlist_out}")
PYEOF
    echo "[$(date)] Joyce TSV parsing complete."

  else
    echo "[$(date)] joyce_annotated_vcf/ found but no recognised annotation files (VCF or TSV)."
    echo "  Contents: $(ls "${JOYCE_DIR}" | tr '\n' ' ')"
    echo "  Falling back to pvar + gnomAD constraint..."
    USE_JOYCE_VCF=0
  fi

else
  echo "[$(date)] joyce_annotated_vcf/ not found — falling back to pvar parsing."
  USE_JOYCE_VCF=0
fi

# ============================================================
# 2. Fallback: parse wes_allchr.pvar + annotate with gnomAD
#    (used only if joyce annotations were not found/parsed)
# ============================================================
if [[ ! -f "${ANNO_OUT}" || ! -s "${ANNO_OUT}" ]]; then
  echo "[$(date)] Falling back: parsing wes_allchr.pvar..."

  # The pvar file has: CHROM POS ID REF ALT [INFO fields]
  # Without functional annotations we can only do approximate classification
  # using variant consequence from the ID field (if rsID) or INFO=CSQ if present.
  # Better alternative: pull gnomAD v4 constraint + LOFTEE from public release.

  PVAR="${UKB_WES_PFILE}.pvar"
  [[ -f "${PVAR}" ]] || { echo "ERROR: wes_allchr.pvar not found: ${PVAR}"; exit 1; }

  echo "  pvar size: $(du -sh ${PVAR})"
  echo "  Peeking at pvar header..."
  head -5 "${PVAR}"

  python - "${PVAR}" "${ANNO_OUT}" "${SETLIST_OUT}" << 'PYEOF'
import sys, gzip, re
from collections import defaultdict

pvar_path   = sys.argv[1]
anno_out    = sys.argv[2]
setlist_out = sys.argv[3]

# pvar format (PLINK2):
# #CHROM  POS  ID  REF  ALT  [FILTER  INFO ...]
# The INFO field may contain CSQ= if this pvar was generated from VEP-annotated VCF.

print(f"  Reading {pvar_path}...")
opener = gzip.open if pvar_path.endswith(".gz") else open

anno_rows   = []
gene_vars   = defaultdict(list)
gene_coords = defaultdict(lambda: {"chrom":"", "min": float("inf"), "max": 0})
n_processed = 0

with opener(pvar_path, "rt") as fh:
    header = None
    for line in fh:
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            header = line.rstrip().lstrip("#").split("\t")
            print(f"  pvar columns: {header}")
            continue
        if header is None:
            continue
        parts = line.rstrip().split("\t")
        chrom  = parts[0].lstrip("chr")
        pos    = parts[1]
        ref    = parts[3]
        alt    = parts[4] if "," not in parts[4] else parts[4].split(",")[0]
        info   = parts[6] if len(parts) > 6 else ""

        # Try to extract CSQ from INFO
        csq_match = re.search(r"CSQ=([^;]+)", info)
        if csq_match:
            csq_str = csq_match.group(1).split(",")[0]
            csq_p   = csq_str.split("|")
            # Basic CSQ field positions (VEP default order):
            # 0=allele, 1=consequence, 2=impact, 3=symbol, 4=gene
            consequence = csq_p[1] if len(csq_p) > 1 else ""
            gene        = csq_p[3] if len(csq_p) > 3 else ""
            lof         = ""
            # Search for LoF field
            for i, val in enumerate(csq_p):
                if val in ("HC", "LC") and i > 5:
                    lof = val
                    break
        else:
            # No CSQ — cannot classify; skip
            continue

        if not gene or gene == ".":
            continue

        if lof == "HC":
            anno = "pLoF"
        elif any(c in consequence for c in
                 ["stop_gained","frameshift","splice_donor","splice_acceptor"]):
            anno = "pLoF"
        elif "missense" in consequence:
            anno = "missense2"  # conservative without REVEL/CADD
        elif "synonymous" in consequence:
            anno = "synonymous"
        else:
            continue

        anno_rows.append(f"{chrom}\t{pos}\t{ref}\t{alt}\t{gene}\t{anno}")
        var_id = f"{chrom}:{pos}:{ref}:{alt}"
        gene_vars[gene].append(var_id)
        p = int(pos)
        gene_coords[gene]["chrom"] = chrom
        gene_coords[gene]["min"]   = min(gene_coords[gene]["min"], p)
        gene_coords[gene]["max"]   = max(gene_coords[gene]["max"], p)
        n_processed += 1
        if n_processed % 500000 == 0:
            print(f"  ... {n_processed:,} variants processed", flush=True)

print(f"  Total: {len(anno_rows):,} variants, {len(gene_vars):,} genes")
with open(anno_out, "w") as fh:
    fh.write("\n".join(anno_rows) + "\n")
with open(setlist_out, "w") as fh:
    for gene, vars_ in sorted(gene_vars.items()):
        c = gene_coords[gene]
        fh.write(f"{gene}\t{c['chrom']}\t{c['min']}\t{c['max']}\t{','.join(vars_)}\n")
print(f"  Written: {anno_out}")
print(f"  Written: {setlist_out}")
PYEOF
  echo "[$(date)] pvar fallback parsing complete."
fi

# ============================================================
# 3. Write mask definition file
# ============================================================
cat > "${MASKS_OUT}" << 'EOF'
mask1	pLoF
mask2	pLoF,missense3
mask3	pLoF,missense3,missense2
EOF

echo ""
echo "[$(date)] Mask definitions written: ${MASKS_OUT}"
cat "${MASKS_OUT}"

# ============================================================
# 4. Validation summary
# ============================================================
echo ""
echo "[$(date)] === Validation ==="
for f in "${ANNO_OUT}" "${SETLIST_OUT}" "${MASKS_OUT}"; do
  if [[ -f "${f}" && -s "${f}" ]]; then
    N=$(wc -l < "${f}")
    echo "  OK  ${f}  (${N} lines)"
  else
    echo "  FAIL  ${f}  (missing or empty)"
  fi
done

if [[ -f "${ANNO_OUT}" ]]; then
  echo ""
  echo "  Annotation category counts:"
  awk '{print $6}' "${ANNO_OUT}" | sort | uniq -c | sort -rn
  echo ""
  echo "  Total genes in set list: $(wc -l < ${SETLIST_OUT})"
fi

echo ""
echo "[$(date)] Burden mask preparation complete."
echo "Ready to submit: sbatch --array=1-22 scripts/agent1_genetics/slurm/03b_burden_test.slurm"
