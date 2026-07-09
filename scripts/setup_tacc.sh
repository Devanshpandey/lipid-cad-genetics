#!/usr/bin/env bash
# Run ONCE on a Lonestar6 LOGIN NODE after cloning the repo.
# Assumes conda envs (regenie_env, r_env, ldsc, pyt_env) already exist.
#
# Uses full binary paths instead of "conda activate" to avoid TACC's
# HOST: unbound variable error in conda activation scripts.
#
# Steps:
#   1. Create project directory structure
#   2. Install system libs (zlib, gmp) into r_env
#   3. Install R genetics packages into r_env (in dependency order)
#   4. Install polars into pyt_env
#   5. Download LDSC LD reference (HapMap3 European LD scores)
#   6. Verify tools

set -uo pipefail   # -e omitted so one failure doesn't abort all steps

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${REPO_DIR}/config/tacc_paths.sh"

# Direct binary paths — avoids conda activate HOST issue on TACC
CONDA_BIN="${CONDA_BASE}/bin/conda"
R_BIN="${CONDA_BASE}/envs/${CONDA_ENV_R}/bin/Rscript"
PIP_PYT="${CONDA_BASE}/envs/${CONDA_ENV_PY}/bin/pip"

echo "============================================"
echo " CAD Genetics — TACC Setup"
echo " Project dir: ${TACC_WORK}"
echo "============================================"

# ---- 1. Directory structure -----------------------------------------------
echo ""
echo "[1/5] Creating directories..."
mkdir -p "${TACC_WORK}"/{data/{reference/{ldsc_eur_w_ld_chr,gtex_v8}},logs}
mkdir -p "${RESULTS_DIR}"/{agent1_genetics/{gwas_summary_stats/{lipids,outcomes},\
phenotypes,regenie_step1,ldsc_munged,ldsc_rg,coloc,mr,finemapping,prs_weights},\
agent2_genes,agent3_networks,agent4_subtypes,final}
mkdir -p "${TACC_SCRATCH}"/{regenie_tmp,ldsc_tmp}
echo "  Done."

# ---- 2. System libs into r_env (zlib → RCurl/XVector; gmp → MendelianRandomization)
echo ""
echo "[2/5] Installing system libraries into r_env via conda..."
"${CONDA_BIN}" install -n "${CONDA_ENV_R}" -c conda-forge zlib gmp -y \
  && echo "  zlib + gmp: OK" \
  || echo "  WARNING: conda install failed — R compilation may still work if libs exist"

# ---- 3. R genetics packages (installed in dependency order) ---------------
echo ""
echo "[3/5] Installing R genetics packages into r_env..."

# Expose conda env's bin so R can find its cross-compilers (x86_64-conda-linux-gnu-cc etc.)
# R's Makeconf hardcodes these; they live in r_env/bin but aren't in PATH unless we add it.
# This avoids `conda activate` (which crashes on TACC with HOST: unbound variable).
CONDA_R_PREFIX="${CONDA_BASE}/envs/${CONDA_ENV_R}"
export PATH="${CONDA_R_PREFIX}/bin:${PATH}"
export CONDA_PREFIX="${CONDA_R_PREFIX}"
export PKG_CONFIG_PATH="${CONDA_R_PREFIX}/lib/pkgconfig:${PKG_CONFIG_PATH:-}"
export LDFLAGS="-L${CONDA_R_PREFIX}/lib ${LDFLAGS:-}"
export CPPFLAGS="-I${CONDA_R_PREFIX}/include ${CPPFLAGS:-}"
export LD_LIBRARY_PATH="${CONDA_R_PREFIX}/lib:${LD_LIBRARY_PATH:-}"

"${R_BIN}" - <<'REOF'
options(repos = c(CRAN = "https://cloud.r-project.org"))

safe_install <- function(pkgs, ...) {
  for (p in pkgs) {
    tryCatch({
      if (!requireNamespace(p, quietly = TRUE)) {
        install.packages(p, quiet = TRUE, ...)
        cat("  installed:", p, "\n")
      } else {
        cat("  already present:", p, "\n")
      }
    }, error = function(e) {
      cat("  FAILED:", p, "—", conditionMessage(e), "\n")
    })
  }
}

safe_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", quiet = TRUE)
  for (p in pkgs) {
    tryCatch({
      if (!requireNamespace(p, quietly = TRUE)) {
        BiocManager::install(p, ask = FALSE, update = FALSE, quiet = TRUE)
        cat("  installed (Bioc):", p, "\n")
      } else {
        cat("  already present:", p, "\n")
      }
    }, error = function(e) {
      cat("  FAILED (Bioc):", p, "—", conditionMessage(e), "\n")
    })
  }
}

# 3a — Core CRAN packages (no special system deps)
cat("--- CRAN core ---\n")
safe_install(c("data.table", "optparse", "ggplot2", "dplyr", "remotes",
               "BiocManager", "coloc", "susieR", "ieugwasr"))

# 3b — Bioconductor base (needs zlib for XVector)
cat("--- Bioconductor base ---\n")
safe_bioc(c("BiocGenerics", "S4Vectors", "IRanges",
            "XVector", "GenomeInfoDb", "GenomicRanges",
            "impute", "preprocessCore", "GO.db", "AnnotationDbi"))

# 3c — WGCNA via BiocManager (install.packages can't find Bioc deps; BiocManager can)
cat("--- WGCNA ---\n")
tryCatch({
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    BiocManager::install("WGCNA", ask = FALSE, update = FALSE, quiet = TRUE)
    cat("  installed: WGCNA\n")
  } else {
    cat("  already present: WGCNA\n")
  }
}, error = function(e) cat("  FAILED: WGCNA —", conditionMessage(e), "\n"))

# 3d — MendelianRandomization (needs gmp)
cat("--- MendelianRandomization ---\n")
safe_install("MendelianRandomization")

# 3e — TwoSampleMR (GitHub, needs ieugwasr)
cat("--- TwoSampleMR ---\n")
tryCatch({
  if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
    remotes::install_github("mrcieu/TwoSampleMR", upgrade = "never", quiet = TRUE)
    cat("  installed: TwoSampleMR\n")
  } else {
    cat("  already present: TwoSampleMR\n")
  }
}, error = function(e) {
  cat("  FAILED: TwoSampleMR —", conditionMessage(e), "\n")
})

# Summary
cat("\n--- Installation summary ---\n")
wanted <- c("data.table","optparse","ggplot2","dplyr","remotes","BiocManager",
            "coloc","susieR","ieugwasr","TwoSampleMR","WGCNA",
            "MendelianRandomization","GenomicRanges","IRanges","XVector")
ok  <- wanted[sapply(wanted, requireNamespace, quietly = TRUE)]
bad <- setdiff(wanted, ok)
cat("OK:", paste(ok, collapse=", "), "\n")
if (length(bad) > 0) cat("MISSING:", paste(bad, collapse=", "), "\n")
REOF

# ---- 4. polars into pyt_env -----------------------------------------------
echo ""
echo "[4/5] Installing polars into pyt_env (direct pip path)..."
"${PIP_PYT}" install polars --quiet \
  && echo "  polars: OK" \
  || echo "  WARNING: polars install failed"

# ---- 5. LDSC LD reference -------------------------------------------------
echo ""
echo "[5/5] Downloading LDSC HapMap3 LD scores (European)..."
REF_DIR="${TACC_WORK}/data/reference"

mkdir -p "${REF_DIR}/ldsc_eur_w_ld_chr"

if [[ ! -f "${REF_DIR}/ldsc_eur_w_ld_chr/22.l2.ldscore.gz" ]]; then
  echo "  Downloading LD scores (~1 GB)..."
  if wget --tries=3 --timeout=120 \
       -O "${REF_DIR}/eur_w_ld_chr.tar.bz2" \
       "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/eur_w_ld_chr.tar.bz2" \
     && tar -xjf "${REF_DIR}/eur_w_ld_chr.tar.bz2" \
          -C "${REF_DIR}/ldsc_eur_w_ld_chr" --strip-components=1; then
    rm -f "${REF_DIR}/eur_w_ld_chr.tar.bz2"
    echo "  LD scores: OK"
  else
    rm -f "${REF_DIR}/eur_w_ld_chr.tar.bz2"
    echo "  LD scores: DOWNLOAD FAILED — retry manually (see script comments)"
  fi
else
  echo "  LD scores: already present"
fi

if [[ ! -f "${REF_DIR}/w_hm3.snplist" ]]; then
  echo "  Downloading HapMap3 SNP list..."
  if wget --tries=3 --timeout=120 \
       -O "${REF_DIR}/w_hm3.snplist.bz2" \
       "https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/w_hm3.snplist.bz2" \
     && bzip2 -d "${REF_DIR}/w_hm3.snplist.bz2"; then
    echo "  HapMap3 SNP list: OK"
  else
    rm -f "${REF_DIR}/w_hm3.snplist.bz2"
    echo "  HapMap3 SNP list: DOWNLOAD FAILED"
  fi
else
  echo "  HapMap3 SNP list: already present"
fi

# ---- Verification ---------------------------------------------------------
echo ""
echo "=== Tool verification ==="
"${PLINK2}" --version 2>&1 | head -1 | sed 's/^/  plink2: /'
"${BCFTOOLS}" --version 2>&1 | head -1 | sed 's/^/  bcftools: /'

# LDSC — use full Python path from ldsc env
LDSC_PY_BIN="${CONDA_BASE}/envs/${CONDA_ENV_LDSC}/bin/python"
"${LDSC_PY_BIN}" "${LDSC_PY}" --help 2>&1 | head -1 | sed 's/^/  ldsc.py: /' || \
  echo "  ldsc.py: check path ${LDSC_PY}"

echo "  UKB GENO step1: $(ls "${UKB_GENO_STEP1}.bed" 2>/dev/null && echo OK || echo MISSING)"
echo "  UKB GENO step2: $(ls "${UKB_GENO_STEP2}.bed" 2>/dev/null && echo OK || echo MISSING)"
echo "  UKB PHENO dir:  $(ls -d "${UKB_PHENO_DIR}" 2>/dev/null && echo OK || echo MISSING)"
echo "  R packages:     $("${R_BIN}" -e 'cat(requireNamespace("TwoSampleMR",quietly=TRUE) && requireNamespace("coloc",quietly=TRUE) && requireNamespace("susieR",quietly=TRUE))' 2>/dev/null)"
echo "  polars:         $("${PIP_PYT}" show polars 2>/dev/null | grep -c Version | sed 's/1/OK/' | sed 's/0/MISSING/')"

cat <<MSG

============================================
 Setup complete (check WARNINGs above).

 Next step — run Agent 1 pipeline:
   cd ${TACC_WORK}
   bash scripts/agent1_genetics/run_all.sh

 Or resume from a specific step:
   bash scripts/agent1_genetics/run_all.sh --start-at 3
============================================
MSG
