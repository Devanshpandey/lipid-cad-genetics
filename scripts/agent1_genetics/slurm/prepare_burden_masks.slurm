#!/bin/bash
#------------------------------------------------------------
# SLURM array: parse joyce_annotated_vcf/ → REGENIE burden annotation files
#
# VCF source (confirmed):
#   annotated_chr_{CHR}_sites_dbNSFP.vcf.gz + .tbi
#   VEP + dbNSFP plugin; NO LOFTEE; NO REVEL/CADD
#
# CSQ FORMAT (verified from chr1 header — exact field indices):
#   0:Allele  1:Consequence  2:IMPACT  3:SYMBOL  7:BIOTYPE  24:CANONICAL
#   43:SIFT(VEP)  44:PolyPhen(VEP)
#   86:LRT_pred   87:LRT_score
#   88:MutationTaster_pred  89:MutationTaster_score
#   90:Polyphen2_HDIV_pred  91:Polyphen2_HDIV_score
#   92:Polyphen2_HVAR_pred  93:Polyphen2_HVAR_score
#   94:SIFT_pred(dbNSFP)    95:SIFT_score
#   96:gnomAD_exomes_AF     98:gnomAD_exomes_POPMAX_AF
#
# Classification (no LOFTEE):
#   pLoF     — IMPACT=HIGH, CANONICAL=YES, protein_coding, POPMAX_AF<0.01
#              (stop_gained, frameshift, splice_donor/acceptor, start/stop_lost)
#   missense3 — IMPACT=MODERATE, CANONICAL, protein_coding, ≥3/4 predictors damaging
#   missense2 — IMPACT=MODERATE, CANONICAL, protein_coding, 2/4 predictors damaging
#   missense1 — IMPACT=MODERATE, CANONICAL, protein_coding, 1/4 predictors damaging
#   synonymous — IMPACT=LOW, synonymous_variant, CANONICAL, protein_coding
#
# Predictors: SIFT_pred=D, PP2_HDIV_pred=D, LRT_pred=D, MutationTaster_pred∈{A,D}
#
# Submit: sbatch --array=1-22 prepare_burden_masks.slurm
# Then:   sbatch --dependency=afterok:<jobid> merge_burden_masks.slurm
#------------------------------------------------------------
#SBATCH -J burden_masks
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 03:00:00
#SBATCH -A OTH26002
#SBATCH -o %x_%A_%a.out
#SBATCH -e %x_%A_%a.err
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=devanshpandey@utexas.edu

set -euo pipefail
# Resolve repo root: SLURM_SUBMIT_DIR when run via sbatch, hardcoded path as fallback
# (BASH_SOURCE[0] points to /var/spool in SLURM jobs — cannot use dirname approach)
_REPO="${SLURM_SUBMIT_DIR:-/work/07880/devansh/lonestar/cad-genetics}"
source "${_REPO}/config/tacc_paths.sh"

CHR="${SLURM_ARRAY_TASK_ID}"

set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_PY}"
set -u

OUT_DIR="${UKB_BURDEN_ANNO_DIR}/by_chr"
mkdir -p "${OUT_DIR}"

VCF="${UKB_WES_JOYCE_VCF}/annotated_chr_${CHR}_sites_dbNSFP.vcf.gz"
ANNO_CHR="${OUT_DIR}/anno_chr${CHR}.txt"
SETLIST_CHR="${OUT_DIR}/setlist_chr${CHR}.txt"

echo "[$(date)] Chr ${CHR} — burden mask annotation"
echo "  Input:  ${VCF}  ($(du -sh ${VCF} | cut -f1))"
echo "  Output: ${ANNO_CHR}"

[[ -f "${VCF}" ]] || { echo "ERROR: VCF not found: ${VCF}"; exit 1; }

if [[ -f "${ANNO_CHR}" && -s "${ANNO_CHR}" ]]; then
  echo "[$(date)] Already done — skipping (delete to rerun)"
  wc -l "${ANNO_CHR}" "${SETLIST_CHR}"
  exit 0
fi

python3 - << PYEOF
import gzip, re, sys
from collections import defaultdict

vcf_path    = "${VCF}"
chrom_label = "${CHR}"
anno_out    = "${ANNO_CHR}"
setlist_out = "${SETLIST_CHR}"

# ── Exact CSQ column indices (verified from chr1 ##INFO=<ID=CSQ header) ──────
# Do NOT change these without re-checking the header.
C_CONSEQUENCE  = 1
C_IMPACT       = 2
C_SYMBOL       = 3
C_BIOTYPE      = 7
C_CANONICAL    = 24
C_LRT_PRED     = 86
C_MTTASTER     = 88    # MutationTaster_pred: A=automatic_disease, D=disease, N=polymorphism, P=poly_auto
C_PP2_HDIV     = 90    # Polyphen2_HDIV_pred: D=probably_damaging, P=possibly, B=benign
C_SIFT_PRED    = 94    # dbNSFP SIFT_pred: D=deleterious, T=tolerated
C_GNOMAD_POPMAX = 98   # gnomAD_exomes_POPMAX_AF

# High-impact consequences that constitute pLoF (no LOFTEE available)
PLOF_CONSEQUENCES = {
    "stop_gained", "frameshift_variant", "splice_donor_variant",
    "splice_acceptor_variant", "start_lost", "stop_lost",
    "transcript_ablation", "splice_donor_5th_base_variant",
    "splice_region_variant",     # conservative inclusion; large effect in exons
}
# Note: splice_region_variant is MODERATE in VEP; include only if it's in an exon
# We flag it as pLoF only when IMPACT=HIGH. For MODERATE splice_region we skip.

MAX_POPMAX_PLOF     = 0.01    # exclude common pLoF (likely annotation artefact)
MAX_POPMAX_MISSENSE = 0.01    # exclude common missense from burden masks

def get_csq(fields, idx):
    return fields[idx] if idx < len(fields) else ""

def is_damaging_sift(v):
    # dbNSFP SIFT_pred: "D" = deleterious; may be multi-value "D;D;T" — any D counts
    return "D" in v.split(";")

def is_damaging_pp2(v):
    # D=probably_damaging counts; P=possibly does not for missense3 but yes for missense2
    return "D" in v.split(";")

def is_possibly_pp2(v):
    return "D" in v.split(";") or "P" in v.split(";")

def is_damaging_lrt(v):
    return "D" in v.split(";")

def is_damaging_mt(v):
    # A=automatic_disease, D=disease_causing both count as damaging
    return any(x in v.split(";") for x in ("A", "D"))

def popmax_af(v):
    if not v or v in (".", ""):
        return 0.0
    try:
        return max(float(x) for x in v.split(";") if x not in (".", ""))
    except ValueError:
        return 0.0

anno_rows    = []
gene_vars    = defaultdict(list)   # gene → [var_ids]
gene_coords  = defaultdict(lambda: {"min": float("inf"), "max": 0})

n_total = n_written = n_skip_af = n_skip_noncanon = 0

with gzip.open(vcf_path, "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        n_total += 1

        cols = line.rstrip().split("\t")
        if len(cols) < 8:
            continue

        chrom  = cols[0].lstrip("chr")
        pos    = cols[1]
        ref    = cols[3]
        alts   = cols[4].split(",")
        info   = cols[7]

        # Extract CSQ= from INFO
        csq_match = re.search(r"CSQ=([^\\t]+)", info)
        if not csq_match:
            continue

        csq_block = csq_match.group(1)

        # Process each ALT allele
        for alt_idx, alt in enumerate(alts):
            if alt in ("*", "<NON_REF>", "."):
                continue

            best_anno = None
            best_gene = None

            # Find canonical transcript entry (CANONICAL=YES, BIOTYPE=protein_coding)
            # VEP outputs one CSQ entry per transcript per allele, separated by comma.
            # The allele in CSQ position 0 may be the VEP-normalised allele (not raw ALT).
            for csq_entry in csq_block.split(","):
                f = csq_entry.split("|")

                # Filter to canonical protein-coding transcripts only
                canonical = get_csq(f, C_CANONICAL)
                biotype   = get_csq(f, C_BIOTYPE)
                if canonical != "YES" or biotype != "protein_coding":
                    n_skip_noncanon += 1
                    continue

                gene        = get_csq(f, C_SYMBOL)
                impact      = get_csq(f, C_IMPACT)
                consequence = get_csq(f, C_CONSEQUENCE)
                af_popmax   = popmax_af(get_csq(f, C_GNOMAD_POPMAX))

                if not gene or gene in (".", ""):
                    continue

                # ── pLoF classification (IMPACT=HIGH, no LOFTEE) ──
                if impact == "HIGH":
                    cons_set = set(c.strip() for c in consequence.split("&"))
                    if cons_set & PLOF_CONSEQUENCES:
                        if af_popmax > MAX_POPMAX_PLOF:
                            n_skip_af += 1
                            continue
                        anno = "pLoF"
                        priority = 5

                # ── Missense classification ──
                elif impact == "MODERATE" and "missense_variant" in consequence:
                    if af_popmax > MAX_POPMAX_MISSENSE:
                        n_skip_af += 1
                        continue

                    sift_pred = get_csq(f, C_SIFT_PRED)
                    pp2_hdiv  = get_csq(f, C_PP2_HDIV)
                    lrt_pred  = get_csq(f, C_LRT_PRED)
                    mt_pred   = get_csq(f, C_MTTASTER)

                    # Count tools with damaging prediction (strict)
                    n_dam = (int(is_damaging_sift(sift_pred)) +
                             int(is_damaging_pp2(pp2_hdiv)) +
                             int(is_damaging_lrt(lrt_pred)) +
                             int(is_damaging_mt(mt_pred)))

                    if n_dam >= 3:
                        anno     = "missense3"
                        priority = 4
                    elif n_dam == 2:
                        anno     = "missense2"
                        priority = 3
                    elif n_dam == 1:
                        anno     = "missense1"
                        priority = 2
                    else:
                        continue  # 0 tools predict damaging → benign

                # ── Synonymous (negative control) ──
                elif impact == "LOW" and "synonymous_variant" in consequence:
                    anno     = "synonymous"
                    priority = 1

                else:
                    continue

                # Keep best annotation (highest priority) across transcripts
                if best_anno is None or priority > best_priority:
                    best_anno     = anno
                    best_gene     = gene
                    best_priority = priority

            if best_anno and best_gene:
                var_id = f"{chrom}:{pos}:{ref}:{alts[0]}"
                anno_rows.append(f"{var_id}\\t{best_gene}\\t{best_anno}")
                gene_vars[best_gene].append(var_id)
                p = int(pos)
                gene_coords[best_gene]["min"] = min(gene_coords[best_gene]["min"], p)
                gene_coords[best_gene]["max"] = max(gene_coords[best_gene]["max"], p)
                n_written += 1

        if n_total % 500_000 == 0:
            print(f"  {n_total:>9,} sites | {n_written:>7,} annotated | "
                  f"{n_skip_af:>6,} AF-filtered | {n_skip_noncanon:>8,} non-canonical",
                  flush=True)

print(f"\n  Chr {chrom_label} complete:")
print(f"    {n_total:,} total variant sites")
print(f"    {n_written:,} annotated variants written")
print(f"    {n_skip_af:,} excluded (AF > threshold)")
print(f"    {n_skip_noncanon:,} skipped (non-canonical / non-coding transcripts)")
print(f"    {len(gene_vars):,} genes")

from collections import Counter
cat_counts = Counter(r.split("\\t")[2] for r in anno_rows)
print("\n  Category breakdown:")
for cat in ["pLoF","missense3","missense2","missense1","synonymous"]:
    print(f"    {cat:12s}: {cat_counts.get(cat, 0):>8,}")

# Write annotation file
with open(anno_out, "w") as fh:
    fh.write("\\n".join(anno_rows) + "\\n")
print(f"\n  Written: {anno_out}")

# Write set list
with open(setlist_out, "w") as fh:
    for gene in sorted(gene_vars):
        c = gene_coords[gene]
        vars_str = ",".join(gene_vars[gene])
        fh.write(f"{gene}\\t{chrom_label}\\t{c['min']}\\t{c['max']}\\t{vars_str}\\n")
print(f"  Written: {setlist_out}")
PYEOF

echo "[$(date)] Chr ${CHR} done."
echo "  Variants: $(wc -l < ${ANNO_CHR})"
echo "  Genes:    $(wc -l < ${SETLIST_CHR})"
