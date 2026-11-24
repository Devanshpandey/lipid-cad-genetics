#!/bin/bash
#=======================================================================
#  run_gwas.sh — one-command GWAS for a single quantitative phenotype
#                UK Biobank on corral / REGENIE / TACC Lonestar6
#=======================================================================
#
#  WHAT IT DOES
#    Give it a phenotype file; it runs the whole common-variant GWAS:
#
#      stage 1   REGENIE step 1   whole-genome null model    1 node
#      stage 2   REGENIE step 2   association testing        22-node array
#      stage 3   merge            sumstats + QC + top hits   1 node
#
#    Measured end to end on UK Biobank LDL-C, N=356,137, 14.8M variants:
#      step 1  2 h 25 min   (with the 987K-SNP shared list; a 164K list
#                            cuts this to well under an hour)
#      step 2  36 min       (slowest chromosome; the array runs in parallel)
#      merge    5 min
#    So roughly 3 h of compute, plus however long LS6 keeps you queued.
#
#    All three are submitted as one SLURM chain with dependencies, so a
#    single command launches everything and you can log off.
#
#  QUICK START
#      ./run_gwas.sh --check                             # verify your setup
#      ./run_gwas.sh --pheno my_trait.txt --covar covariates.txt
#
#    my_trait.txt is whitespace- or tab-delimited WITH a header row:
#
#      FID      IID      MYTRAIT
#      <ID>  <ID>  4.21
#      <ID>  <ID>  NA
#
#    Missing values are NA. Samples not listed are simply dropped.
#    If the file has exactly three columns, --trait can be omitted.
#
#  WHAT YOU GET   (in ./<TRAIT>/ unless you pass --out)
#      <TRAIT>.gwas.tsv.gz    genome-wide summary statistics
#      <TRAIT>.tophits.tsv    independent genome-wide-significant loci
#      <TRAIT>.summary.txt    N, lambda_GC, hit count, full provenance
#      run_config.sh          every resolved setting for this run
#      logs/                  all SLURM stdout/stderr for this run
#
#  WHAT IS FIXED vs WHAT IS YOURS
#    Fixed:  the shared UKB genotypes and QC'd sample list on corral.
#            Nobody needs to copy or re-derive these.
#    Yours:  allocation, scratch, output directory and REGENIE install are
#            all detected from your own environment — no paths belonging to
#            one person are baked in. Phenotype and covariate files are
#            per-project, so you pass them per run.
#
#  RESUMING
#    Re-running the same command skips work already finished (step-1
#    predictions, per-chromosome results). Delete the output directory to
#    force a clean run. Rerun a single stage with --stage (see --help).
#
#  QUANTITATIVE TRAITS ONLY. A binary phenotype is detected and refused:
#  binary GWAS needs REGENIE --bt with Firth correction, not this script.
#=======================================================================

set -euo pipefail
GWAS_VERSION="1.2"

#=======================================================================
#  SHARED LAB DATA ON CORRAL — same for everyone, read-only
#=======================================================================
UKB_GENO_ROOT="/path/to/ukb_genotypes"

# Step 1 (null model) uses the MAF>0.01 merged file; step 2 (association)
# uses the MAF>0.001 merged file. Override per run with --geno-step1/2.
GENO_STEP1="${UKB_GENO_ROOT}/merged_maf0.01_biallel_bbf_400k/merged_sub_chrom_maf0.01"
GENO_STEP2="${UKB_GENO_ROOT}/merged_maf0.001_biallel_bbf_400k/merged_sub_chrom_maf0.001"

# QC-passed white-British sample list, already in REGENIE FID/IID format.
KEEP_DEFAULT="${UKB_GENO_ROOT}/geno_qc_eids_400k_white_british_2023-11-29.txt"

# Lab-wide LD-pruned SNP list for the step-1 null model (~987K SNPs, 99% of
# which are in the step-1 genotypes). Using it means nobody needs PLINK2.
# Override with --pruned-snps; a smaller set makes step 1 faster.
PRUNED_SNPS_SHARED="/path/to/ukb/GWAS_INPUT_FILES/prune_list_maf0.01.prune.in"

SAMPLE_FILE=""        # BGEN .sample file, only if you point at BGEN input

#=======================================================================
#  YOUR ENVIRONMENT — all auto-detected; override with flags if needed
#=======================================================================
ACCOUNT=""            # "" = first allocation from sacctmgr; --account to pick
PARTITION="normal"    # LS6 CPU queue
EMAIL=""              # "" = no mail; --email to enable
CONDA_ENV="regenie_env"   # only consulted if regenie is not already on PATH
REGENIE_BIN=""        # --regenie to point at a specific binary
PLINK2_BIN=""         # --plink2  to point at a specific binary
OUT_ROOT=""           # "" = current directory
SCRATCH_ROOT=""       # "" = $SCRATCH
PRUNED_SNPS=""        # "" = cached under $SCRATCH, built once with PLINK2

#=======================================================================
#  ANALYSIS DEFAULTS
#=======================================================================
MIN_MAC=20            # minimum minor-allele count (common-variant GWAS)
MAF=""                # optional MAF floor, applied when merging
APPLY_RINT=1          # rank inverse-normal transform in step 2
BSIZE_STEP1=1000
BSIZE_STEP2=400
GW_SIG="5e-8"         # genome-wide significance threshold
CLUMP_KB=1000         # distance pruning window for top hits (kb)
CHROMS="1 22"         # first and last chromosome to analyse

TIME_STEP1="06:00:00"
TIME_STEP2="06:00:00"
TIME_MERGE="02:00:00"
NTASK_STEP1=64
NTASK_STEP2=16
NTASK_MERGE=8

#=======================================================================
#  HELP
#=======================================================================
usage() {
cat <<'__HELPTEXT__'
run_gwas.sh — single-phenotype quantitative GWAS, UK Biobank on TACC

SYNOPSIS
  ./run_gwas.sh --check                       verify your environment
  ./run_gwas.sh --pheno FILE [--covar FILE] [--trait NAME] [options]

REQUIRED
  --pheno FILE        Phenotype file: header row, columns FID IID <TRAIT> [...]
                      Whitespace- or tab-delimited. Missing = NA.

STRONGLY RECOMMENDED
  --covar FILE        Covariate file (FID IID <covariates>), e.g. age, sex,
                      PC1-PC10. Without it nothing corrects for population
                      structure and your results will be inflated.

COMMON OPTIONS
  --trait NAME        Phenotype column to analyse.
                      Optional if the file has exactly 3 columns.
  --out DIR           Output directory.  Default: ./<TRAIT>
  --covar-cols LIST   Use only these covariate columns, comma-separated.
                      REGENIE ranges work:  age,sex,PC{1:10}
  --cat-covar LIST    Covariates to treat as categorical, comma-separated.
  --keep FILE         Restrict to these samples (FID IID, no header).
                      Default: the QC'd white-British list on corral.
                      Pass "" to use every genotyped sample.

DATA OPTIONS  (defaults point at the shared corral genotypes)
  --geno-step1 PREFIX   Genotypes for the null model. Default: corral MAF>0.01
  --geno-step2 PREFIX   Genotypes to test.           Default: corral MAF>0.001
                        bed / pgen / bgen auto-detected.
  --pruned-snps FILE    LD-pruned SNP list for step 1. Defaults to the shared
                        ~987K list on corral, so nothing has to be built.
                        A smaller set makes step 1 faster; pass your own here.
                        If the shared list is unreachable, one is built with
                        PLINK2 and cached under $SCRATCH.
  --chroms "A B"        First and last chromosome.   Default: "1 22"

ANALYSIS OPTIONS
  --min-mac INT       Minimum minor-allele count.    Default: 20
  --maf FLOAT         MAF floor applied when merging (e.g. 0.01). Variants
                      below it are dropped from the final sumstats; the
                      per-chromosome files keep everything.
  --no-rint           Do NOT rank-inverse-normal transform the phenotype.
                      Default is to apply it — right for skewed biomarkers;
                      turn off only if the trait is already transformed.
  --force-binary      Run even if the phenotype looks binary. Gives a linear
                      approximation, normally the wrong analysis.

ENVIRONMENT OPTIONS  (normally auto-detected)
  --account CODE      SLURM allocation. Default: your first from sacctmgr.
  --partition NAME    SLURM partition.  Default: normal
  --email ADDR        Send job mail. Default: none.
  --threads N         Threads per job. Default: all cores the job holds.
                      Measured on chr22: 16 threads 261s, 128 threads 259s,
                      byte-identical results. Step 2 is I/O-bound, so extra
                      threads buy little — but do not let it fall to 1.
  --regenie PATH      REGENIE binary.   Default: PATH, else conda regenie_env
  --plink2 PATH       PLINK2 binary.    Only needed to build a pruned SNP list
  --time-step1 HH:MM:SS   --time-step2 HH:MM:SS   --time-merge HH:MM:SS

CONTROL
  --stage STAGE       all (default) | step1 | step2 | merge
                      Rerun one stage of an existing run:
                        ./run_gwas.sh --stage merge --config <out>/run_config.sh
                      For step2 outside an array job, set CHR=N first.
  --config FILE       Reuse a previous run's run_config.sh.
  --dry-run           Validate everything and print the sbatch commands only.
  --version           Print version.
  -h, --help          This message.

EXAMPLES
  # Typical run — 3-column phenotype file, trait name detected
  ./run_gwas.sh --pheno pheno_tghdl.txt --covar covariates.txt

  # Named column of a multi-trait file, subset of covariates
  ./run_gwas.sh --pheno pheno_lipids.txt --trait LDL_C \
                --covar covariates.txt --covar-cols age,sex,PC{1:10}

  # Pre-residualised trait, no transform, output on scratch
  ./run_gwas.sh --pheno resid.txt --trait LVEF_resid --no-rint \
                --out $SCRATCH/gwas/lvef

  # Test the exome data instead of the array data
  ./run_gwas.sh --pheno t.txt --covar c.txt \
     --geno-step2 /path/to/ukb_wes/wes_allchr

  # Check inputs without submitting anything
  ./run_gwas.sh --pheno my_trait.txt --covar covariates.txt --dry-run

TYPICAL RUNTIME   ~1.5 h wall clock end to end at N ~ 400,000.
__HELPTEXT__
}

#=======================================================================
#  UTILITIES
#=======================================================================
say()  { printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
die()  { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
warn() { printf 'WARNING: %s\n' "$*" >&2; }

need_file() {
  [[ -f "$1" ]] || die "$2 not found: $1"
  [[ -r "$1" ]] || die "$2 is not readable: $1"
}

# Absolute path for a file that is already known to exist.
abspath() { printf '%s/%s\n' "$(cd "$(dirname "$1")" && pwd)" "$(basename "$1")"; }

#=======================================================================
#  ENVIRONMENT DETECTION — nothing here is specific to any one user
#=======================================================================
autodetect() {
  [[ -z "$SCRATCH_ROOT" ]] && SCRATCH_ROOT="${SCRATCH:-${TMPDIR:-/tmp}}"
  [[ -z "$OUT_ROOT"     ]] && OUT_ROOT="$PWD"
  if [[ -z "$PRUNED_SNPS" ]]; then
    if [[ -f "$PRUNED_SNPS_SHARED" ]]; then
      PRUNED_SNPS="$PRUNED_SNPS_SHARED"
    else
      PRUNED_SNPS="${SCRATCH_ROOT}/gwas_cache/$(basename "$GENO_STEP1")_pruned.snplist"
    fi
  fi
  if [[ -z "$ACCOUNT" ]] && command -v sacctmgr >/dev/null 2>&1; then
    ACCOUNT=$(sacctmgr show assoc user="${USER}" format=account -nP 2>/dev/null \
              | sort -u | grep -v '^default$' | head -1 || true)
  fi
  return 0
}

# Make REGENIE and PLINK2 callable. Safe to call more than once.
# Prefers whatever is already on PATH; falls back to a conda env.
load_tools() {
  export HOST="${HOST:-${HOSTNAME:-$(hostname)}}"   # conda 4.x refs $HOST under set -u
  REGENIE_CMD="${REGENIE_BIN:-$(command -v regenie 2>/dev/null || true)}"
  PLINK2_CMD="${PLINK2_BIN:-$(command -v plink2 2>/dev/null || true)}"

  if [[ -z "$REGENIE_CMD" && -n "$CONDA_ENV" && -n "${CONDA_EXE:-}" ]]; then
    local base; base="$(dirname "$(dirname "$CONDA_EXE")")"
    if [[ -f "${base}/etc/profile.d/conda.sh" ]]; then
      set +u
      # shellcheck disable=SC1091
      source "${base}/etc/profile.d/conda.sh"
      conda activate "${CONDA_ENV}" 2>/dev/null || true
      set -u
      REGENIE_CMD="$(command -v regenie 2>/dev/null || true)"
      PLINK2_CMD="${PLINK2_CMD:-$(command -v plink2 2>/dev/null || true)}"
    fi
  fi
}

# How many threads REGENIE should use.
#
# NOT `nproc`. TACC sets OMP_NUM_THREADS=1 site-wide and nproc honours it,
# so nproc reports 1 on both login and compute nodes even though the job
# owns the whole 128-core node (measured: affinity list 0-127). Using it
# made a chr22 step-2 run take 619s instead of 231s.
#
# NOT `nproc --all` either: that reports the physical machine, which is
# 256 on an LS6 login node and 128 on a compute node.
#
# SLURM_CPUS_ON_NODE is what the job actually owns here, so prefer it.
n_threads() {
  local n="${THREADS:-}"                     # --threads wins if given
  [[ -z "$n" ]] && n="${SLURM_CPUS_ON_NODE:-}"
  [[ -z "$n" ]] && n="${SLURM_NTASKS:-}"
  [[ -z "$n" ]] && n="$(nproc --all 2>/dev/null || echo 16)"
  printf '%s\n' "$n"
}

# Genotype arguments for REGENIE, format auto-detected from what exists.
GENO_ARGS=()
set_geno_args() {
  local p="$1" what="$2"
  GENO_ARGS=()
  if   [[ -f "${p}.bed"  ]]; then GENO_ARGS=(--bed  "$p")
  elif [[ -f "${p}.pgen" ]]; then GENO_ARGS=(--pgen "$p")
  elif [[ -f "${p}.bgen" ]]; then
    GENO_ARGS=(--bgen "${p}.bgen")
    [[ -n "$SAMPLE_FILE" ]] && GENO_ARGS+=(--sample "$SAMPLE_FILE")
  else
    die "no genotype files for ${what}: ${p}
     Expected one of ${p}.bed / ${p}.pgen / ${p}.bgen"
  fi
}

geno_format() {   # for reporting only
  local p="$1"
  if   [[ -f "${p}.bed"  ]]; then echo "PLINK bed"
  elif [[ -f "${p}.pgen" ]]; then echo "PLINK2 pgen"
  elif [[ -f "${p}.bgen" ]]; then echo "BGEN"
  else echo "MISSING"; fi
}

#=======================================================================
#  --check — verify the setup before anyone waits on a queue
#=======================================================================
run_check() {
  autodetect
  load_tools
  local fail=0
  local ok="  OK   " no="  MISS "

  echo
  echo "  run_gwas.sh ${GWAS_VERSION} — setup check for ${USER}"
  echo "  ---------------------------------------------------------------"

  if [[ -n "$REGENIE_CMD" ]]; then
    # `|| true` guards the pipeline: a multi-line version banner would
    # otherwise SIGPIPE regenie, and pipefail would abort the check.
    local ver; ver="$("$REGENIE_CMD" --version 2>&1 | head -1 | tr -d '\n' || true)"
    printf '%s regenie          %s (%s)\n' "$ok" "$REGENIE_CMD" "${ver:-version unknown}"
  else
    printf '%s regenie          not found\n' "$no"
    echo "         Activate an env that has it, or pass --regenie /path/to/regenie."
    echo "         This script also tries the conda env '${CONDA_ENV}' automatically."
    fail=1
  fi

  if command -v sbatch >/dev/null 2>&1; then
    printf '%s sbatch           %s\n' "$ok" "$(command -v sbatch)"
  else
    printf '%s sbatch           not found — are you on an LS6 login node?\n' "$no"; fail=1
  fi

  echo
  local f
  for f in GENO_STEP1 GENO_STEP2; do
    local prefix="${!f}" fmt
    fmt="$(geno_format "$prefix")"
    if [[ "$fmt" == "MISSING" ]]; then
      printf '%s %-16s no .bed/.pgen/.bgen at %s\n' "$no" "$f" "$prefix"
      echo "         Check that you have read access to corral."
      fail=1
    else
      printf '%s %-16s %s  [%s]\n' "$ok" "$f" "$(basename "$prefix")" "$fmt"
    fi
  done

  if [[ -f "$KEEP_DEFAULT" ]]; then
    printf '%s keep list        corral white-British (%s samples)\n' "$ok" "$(wc -l < "$KEEP_DEFAULT")"
  else
    printf '%s keep list        not readable: %s\n' "$no" "$KEEP_DEFAULT"; fail=1
  fi

  echo
  if [[ "$PRUNED_SNPS" == "$PRUNED_SNPS_SHARED" && -f "$PRUNED_SNPS" ]]; then
    printf '%s pruned SNPs      shared corral list (%s SNPs)\n' "$ok" "$(wc -l < "$PRUNED_SNPS")"
  elif [[ -f "$PRUNED_SNPS" ]]; then
    printf '%s pruned SNPs      %s (%s SNPs)\n' "$ok" "$PRUNED_SNPS" "$(wc -l < "$PRUNED_SNPS")"
  elif [[ -n "$PLINK2_CMD" ]]; then
    printf '%s pruned SNPs      will be built once with %s\n' "$ok" "$PLINK2_CMD"
    printf '                    and cached at %s\n' "$PRUNED_SNPS"
  else
    printf '%s pruned SNPs      no cached list, and PLINK2 not found\n' "$no"
    echo "         Pass --plink2 /path/to/plink2 to build one (one-off), or"
    echo "         --pruned-snps to reuse a list a labmate already built."
    fail=1
  fi

  echo
  printf '  --   account          %s\n' "${ACCOUNT:-<none detected — pass --account>}"
  printf '  --   partition        %s\n' "$PARTITION"
  if [[ -d "$SCRATCH_ROOT" && -w "$SCRATCH_ROOT" ]]; then
    printf '%s scratch          %s\n' "$ok" "$SCRATCH_ROOT"
  else
    printf '%s scratch          not writable: %s\n' "$no" "$SCRATCH_ROOT"; fail=1
  fi
  printf '  --   output root      %s\n' "$OUT_ROOT"
  echo
  echo "  Covariates are per-project — pass --covar with age, sex and PCs."
  echo "  ---------------------------------------------------------------"
  if (( fail )); then
    echo "  Not ready. Fix the MISS lines above, then re-run --check."
    echo
    return 1
  fi
  echo "  Ready. Try:"
  echo "    ./run_gwas.sh --pheno <file> --covar <file> --dry-run"
  echo
}

#=======================================================================
#  ARGUMENT PARSING
#=======================================================================
PHENO=""; TRAIT=""; OUTDIR=""; COVAR=""; COVAR_COLS=""; CAT_COVAR=""
KEEP=""; STAGE="all"; CONFIG=""; DRY_RUN=0; FORCE_BINARY=0; MODE="run"; THREADS=""; PARTITION_SET=0
REGENIE_CMD=""; PLINK2_CMD=""
# Track an explicit --keep so that passing "" can switch OFF the corral
# default. ${VAR:-default} cannot express that: it fires on empty too.
KEEP_SET=0

# Resolve our own path BEFORE sourcing a config, so a saved config can
# override it. Inside a SLURM job BASH_SOURCE points at the spool copy,
# which is why the real path is persisted into run_config.sh.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --check)        MODE="check";      shift   ;;
    --pheno)        PHENO="$2";        shift 2 ;;
    --trait)        TRAIT="$2";        shift 2 ;;
    --out)          OUTDIR="$2";       shift 2 ;;
    --covar)        COVAR="$2";        shift 2 ;;
    --covar-cols)   COVAR_COLS="$2";   shift 2 ;;
    --cat-covar)    CAT_COVAR="$2";    shift 2 ;;
    --keep)         KEEP="$2"; KEEP_SET=1; shift 2 ;;
    --geno-step1)   GENO_STEP1="$2";   shift 2 ;;
    --geno-step2)   GENO_STEP2="$2";   shift 2 ;;
    --pruned-snps)  PRUNED_SNPS="$2";  shift 2 ;;
    --chroms)       CHROMS="$2";       shift 2 ;;
    --maf)          MAF="$2";          shift 2 ;;
    --min-mac)      MIN_MAC="$2";      shift 2 ;;
    --no-rint)      APPLY_RINT=0;      shift   ;;
    --force-binary) FORCE_BINARY=1;    shift   ;;
    --account)      ACCOUNT="$2";      shift 2 ;;
    --partition)    PARTITION="$2"; PARTITION_SET=1; shift 2 ;;
    --email)        EMAIL="$2";        shift 2 ;;
    --regenie)      REGENIE_BIN="$2";  shift 2 ;;
    --plink2)       PLINK2_BIN="$2";   shift 2 ;;
    --threads)      THREADS="$2";      shift 2 ;;
    --time-step1)   TIME_STEP1="$2";   shift 2 ;;
    --time-step2)   TIME_STEP2="$2";   shift 2 ;;
    --time-merge)   TIME_MERGE="$2";   shift 2 ;;
    --stage)        STAGE="$2";        shift 2 ;;
    --config)       CONFIG="$2";       shift 2 ;;
    --dry-run)      DRY_RUN=1;         shift   ;;
    --version)      echo "run_gwas.sh ${GWAS_VERSION}"; exit 0 ;;
    -h|--help)      usage; exit 0 ;;
    *) die "Unknown option: $1  (run with --help)" ;;
  esac
done

# A run_config.sh from a previous submission restores every setting.
# Flags given on this command line still win over it, so a stage can be
# rerun with, say, a different --threads or --account without editing the
# frozen config. Stash the CLI values, source, then put them back.
if [[ -n "$CONFIG" ]]; then
  need_file "$CONFIG" "Config file"
  _cli_threads="$THREADS"; _cli_account="$ACCOUNT"; _cli_partition="$PARTITION"
  _cli_regenie="$REGENIE_BIN"; _cli_plink2="$PLINK2_BIN"
  # shellcheck disable=SC1090
  source "$CONFIG"
  [[ -n "$_cli_threads" ]] && THREADS="$_cli_threads"
  [[ -n "$_cli_account" ]] && ACCOUNT="$_cli_account"
  [[ -n "$_cli_regenie" ]] && REGENIE_BIN="$_cli_regenie"
  [[ -n "$_cli_plink2"  ]] && PLINK2_BIN="$_cli_plink2"
  (( PARTITION_SET )) && PARTITION="$_cli_partition"
  KEEP_SET=1     # a saved config already holds the final value
fi

CHR_FIRST="${CHROMS%% *}"
CHR_LAST="${CHROMS##* }"

#=======================================================================
#  STAGE: step1 — REGENIE null model (whole-genome regression)
#=======================================================================
build_pruned_snps() {
  # One-off: LD-prune the step-1 genotypes down to the ~100-500K
  # independent common SNPs REGENIE wants for the null model. Cached
  # under scratch, so this cost is paid once per dataset, not per trait.
  [[ -n "$PLINK2_CMD" ]] || die "no cached LD-pruned SNP list and PLINK2 not available.
     Pass --plink2 /path/to/plink2, or --pruned-snps to reuse an existing list."
  local threads; threads="$(n_threads)"
  local tmp="${PRUNED_SNPS%.*}.tmp"
  mkdir -p "$(dirname "$PRUNED_SNPS")"
  say "building LD-pruned SNP list (one-off) -> ${PRUNED_SNPS}"
  set_geno_args "$GENO_STEP1" "step 1"
  local -a keep_arg=()
  [[ -n "$KEEP" ]] && keep_arg=(--keep "$KEEP")
  "$PLINK2_CMD" "${GENO_ARGS[@]}" "${keep_arg[@]+"${keep_arg[@]}"}" \
    --maf 0.01 --geno 0.01 --hwe 1e-15 \
    --indep-pairwise 1000 100 0.1 \
    --threads "$threads" --out "$tmp"
  mv "${tmp}.prune.in" "$PRUNED_SNPS"
  rm -f "${tmp}".*
  say "pruned SNP list: $(wc -l < "$PRUNED_SNPS") SNPs"
}

run_step1() {
  autodetect; load_tools
  [[ -n "$REGENIE_CMD" ]] || die "regenie not found — run ./run_gwas.sh --check"
  local threads; threads="$(n_threads)"
  export OMP_NUM_THREADS="$threads"   # TACC defaults this to 1 site-wide
  local tmpdir="${SCRATCH_ROOT}/regenie_tmp/${TRAIT}_${SLURM_JOB_ID:-$$}"
  mkdir -p "$tmpdir" "${OUTDIR}/step1"

  if [[ -f "${OUTDIR}/step1/step1_${TRAIT}_pred.list" ]]; then
    say "step 1 predictions already exist — skipping"
    return 0
  fi

  [[ -f "$PRUNED_SNPS" ]] || build_pruned_snps

  say "REGENIE step 1 — null model for ${TRAIT} (${threads} threads)"
  set_geno_args "$GENO_STEP1" "step 1"

  local -a args=(
    --step 1
    "${GENO_ARGS[@]}"
    --extract    "${PRUNED_SNPS}"
    --phenoFile  "${PHENO}"
    --phenoCol   "${TRAIT}"
    --bsize      "${BSIZE_STEP1}"
    --lowmem
    --lowmem-prefix "${tmpdir}/step1"
    --threads    "${threads}"
    --gz
    --out        "${OUTDIR}/step1/step1_${TRAIT}"
  )
  [[ -n "$KEEP"       ]] && args+=(--keep "$KEEP")
  [[ -n "$COVAR"      ]] && args+=(--covarFile "$COVAR")
  [[ -n "$COVAR_COLS" ]] && args+=(--covarColList "$COVAR_COLS")
  [[ -n "$CAT_COVAR"  ]] && args+=(--catCovarList "$CAT_COVAR")

  "$REGENIE_CMD" "${args[@]}"
  rm -rf "$tmpdir"

  [[ -f "${OUTDIR}/step1/step1_${TRAIT}_pred.list" ]] \
    || die "step 1 finished but wrote no prediction list"
  say "step 1 complete -> ${OUTDIR}/step1/step1_${TRAIT}_pred.list"
}

#=======================================================================
#  STAGE: step2 — association test, one chromosome per array task
#=======================================================================
run_step2() {
  autodetect; load_tools
  [[ -n "$REGENIE_CMD" ]] || die "regenie not found — run ./run_gwas.sh --check"
  local chr="${SLURM_ARRAY_TASK_ID:-${CHR:-}}"
  [[ -n "$chr" ]] || die "step2 needs SLURM_ARRAY_TASK_ID (or CHR=N) set"
  local threads; threads="$(n_threads)"
  export OMP_NUM_THREADS="$threads"   # TACC defaults this to 1 site-wide
  mkdir -p "${OUTDIR}/per_chr"

  local out="${OUTDIR}/per_chr/chr${chr}"
  if [[ -f "${out}_${TRAIT}.regenie.gz" ]]; then
    say "chr${chr} already done — skipping"
    return 0
  fi

  need_file "${OUTDIR}/step1/step1_${TRAIT}_pred.list" "Step-1 prediction list"
  say "REGENIE step 2 — ${TRAIT} chr${chr} (${threads} threads)"
  set_geno_args "$GENO_STEP2" "step 2"

  local -a args=(
    --step 2
    "${GENO_ARGS[@]}"
    --phenoFile  "${PHENO}"
    --phenoCol   "${TRAIT}"
    --pred       "${OUTDIR}/step1/step1_${TRAIT}_pred.list"
    --chr        "${chr}"
    --bsize      "${BSIZE_STEP2}"
    --minMAC     "${MIN_MAC}"
    --threads    "${threads}"
    --gz
    --out        "${out}"
  )
  [[ "$APPLY_RINT" -eq 1 ]] && args+=(--apply-rint)
  [[ -n "$KEEP"          ]] && args+=(--keep "$KEEP")
  [[ -n "$COVAR"         ]] && args+=(--covarFile "$COVAR")
  [[ -n "$COVAR_COLS"    ]] && args+=(--covarColList "$COVAR_COLS")
  [[ -n "$CAT_COVAR"     ]] && args+=(--catCovarList "$CAT_COVAR")

  "$REGENIE_CMD" "${args[@]}"
  say "chr${chr} complete -> ${out}_${TRAIT}.regenie.gz"
}

#=======================================================================
#  STAGE: merge — concatenate chromosomes, add P, QC, extract top hits
#=======================================================================
run_merge() {
  autodetect
  local merged="${OUTDIR}/${TRAIT}.gwas.tsv"
  local tophits="${OUTDIR}/${TRAIT}.tophits.tsv"
  local summary="${OUTDIR}/${TRAIT}.summary.txt"

  say "Merging chromosomes ${CHR_FIRST}-${CHR_LAST} for ${TRAIT}"

  #--- confirm every chromosome landed ---------------------------------
  local missing=() chr
  for chr in $(seq "$CHR_FIRST" "$CHR_LAST"); do
    [[ -f "${OUTDIR}/per_chr/chr${chr}_${TRAIT}.regenie.gz" ]] || missing+=("$chr")
  done
  if (( ${#missing[@]} > 0 )); then
    die "missing chromosome results: ${missing[*]}
     Rerun just those chromosomes with:
       sbatch --array=$(IFS=,; echo "${missing[*]}") ${ACCOUNT:+-A ${ACCOUNT} }-p ${PARTITION} \\
         -N 1 -n ${NTASK_STEP2} -t ${TIME_STEP2} -J ${TRAIT}_s2 \\
         -o ${OUTDIR}/logs/step2_%A_%a.out -e ${OUTDIR}/logs/step2_%A_%a.err \\
         ${SCRIPT_PATH} --stage step2 --config ${OUTDIR}/run_config.sh"
  fi

  #--- concatenate ------------------------------------------------------
  # REGENIE step-2 columns (space-delimited):
  #   1 CHROM  2 GENPOS  3 ID  4 ALLELE0  5 ALLELE1  6 A1FREQ  7 N
  #   8 TEST   9 BETA   10 SE  11 CHISQ  12 LOG10P  13 EXTRA
  # ALLELE1 is the tested/effect allele; BETA and A1FREQ refer to it.
  printf 'CHR\tBP\tSNP\tEA\tNEA\tEAF\tN\tBETA\tSE\tCHISQ\tLOG10P\tP\n' > "$merged"
  for chr in $(seq "$CHR_FIRST" "$CHR_LAST"); do
    gunzip -c "${OUTDIR}/per_chr/chr${chr}_${TRAIT}.regenie.gz" \
      | awk -v maf="${MAF}" 'NR>1 && $12!="NA" {
          if (maf != "") { f = $6 + 0; if (f > 0.5) f = 1 - f; if (f < maf + 0) next }
          p = 10 ^ (-$12)
          printf "%s\t%s\t%s\t%s\t%s\t%.6g\t%d\t%.6g\t%.6g\t%.6g\t%.6g\t%.6g\n",
                 $1, $2, $3, $5, $4, $6, $7, $9, $10, $11, $12, p
        }' >> "$merged"
  done

  local nsnp; nsnp=$(( $(wc -l < "$merged") - 1 ))
  (( nsnp > 0 )) || die "merged file is empty — check ${OUTDIR}/logs/step2_*.err"

  #--- genomic inflation factor ----------------------------------------
  # lambda_GC = median(chi-square) / 0.4549364  (median of the 1-df null)
  # No early awk exit: that would SIGPIPE sort and trip pipefail.
  local lambda
  lambda=$(tail -n +2 "$merged" | cut -f10 | sort -g -S 2G \
           | awk -v n="$nsnp" 'NR == int((n+1)/2) { m = $1 }
                               END { printf "%.4f", m / 0.4549364 }')

  #--- max N and significant-hit count (single pass, no pipes) ---------
  local maxn nsig
  maxn=$(awk -F'\t' 'NR>1 && $7+0 > m { m = $7+0 } END { print m+0 }' "$merged")
  nsig=$(awk -F'\t' -v thr="$GW_SIG" 'NR>1 && $12+0 < thr+0 { c++ } END { print c+0 }' "$merged")

  #--- independent top hits (greedy distance pruning) -------------------
  head -1 "$merged" > "$tophits"
  awk -F'\t' -v thr="$GW_SIG" 'NR>1 && $12+0 < thr+0' "$merged" \
    | sort -t$'\t' -k11,11gr \
    | awk -F'\t' -v w=$(( CLUMP_KB * 1000 )) '
        {
          keep = 1
          for (i = 1; i <= n; i++)
            if (c[i] == $1 && ($2 - p[i] < w && p[i] - $2 < w)) { keep = 0; break }
          if (keep) { n++; c[n] = $1; p[n] = $2; print }
        }' \
    | sort -t$'\t' -k1,1n -k2,2n >> "$tophits"
  local nloci; nloci=$(( $(wc -l < "$tophits") - 1 ))

  #--- summary ----------------------------------------------------------
  {
    echo "GWAS summary — ${TRAIT}"
    echo "generated          $(date '+%Y-%m-%d %H:%M:%S') by ${USER:-unknown}"
    echo "run_gwas.sh        ${GWAS_VERSION}"
    echo
    echo "phenotype file     ${PHENO}"
    echo "phenotype column   ${TRAIT}"
    echo "covariate file     ${COVAR:-<none>}"
    echo "covariate columns  ${COVAR_COLS:-<all columns in file>}"
    echo "categorical covar  ${CAT_COVAR:-none}"
    echo "sample keep list   ${KEEP:-<all genotyped samples>}"
    echo "genotypes          ${GENO_STEP2}"
    echo "null-model geno    ${GENO_STEP1}"
    echo "pruned SNP list    ${PRUNED_SNPS}"
    echo "rank-inv-normal    $( ((APPLY_RINT)) && echo yes || echo no )"
    echo "min MAC            ${MIN_MAC}"
    echo "MAF floor          ${MAF:-none}"
    echo "chromosomes        ${CHR_FIRST}-${CHR_LAST}"
    echo
    echo "max N tested       ${maxn}"
    echo "variants tested    ${nsnp}"
    echo "lambda_GC          ${lambda}"
    echo "P < ${GW_SIG}        ${nsig} variants"
    echo "independent loci   ${nloci}  (${CLUMP_KB} kb distance pruning)"
    echo
    echo "NOTES"
    echo "  EA is the effect allele (REGENIE ALLELE1); BETA and EAF are"
    echo "  relative to EA. P underflows to 0 below ~1e-308 — use LOG10P"
    echo "  for the strongest variants. lambda_GC well above 1 with a"
    echo "  polygenic trait is expected; check it against the LDSC"
    echo "  intercept before calling it confounding."
  } > "$summary"

  gzip -f "$merged"

  say "merge complete"
  echo
  cat "$summary"
  echo
  say "sumstats : ${merged}.gz"
  say "top hits : ${tophits}"
}

# Submit a job and echo just its numeric id.
#
# TACC's sbatch is a wrapper that prints a multi-line banner to STDOUT
# ("Welcome to Lonestar6", "Verifying valid submit host...OK", ...), so
# `sbatch ... | awk '{print $NF}'` returns the last field of EVERY line,
# not the job id — which then gets passed to --dependency as junk and the
# chain fails. Match the canonical SLURM line, and fail loudly on garbage.
submit_job() {
  local out jid
  out="$(sbatch "$@" 2>&1)" || die "sbatch failed:
${out}"
  jid="$(printf '%s\n' "$out" | awk '/Submitted batch job/ { print $NF }' | tail -1)"
  [[ "$jid" =~ ^[0-9]+$ ]] || die "could not read a job id from sbatch output:
${out}"
  printf '%s\n' "$jid"
}

#=======================================================================
#  STAGE: all — validate inputs, then submit the SLURM chain
#=======================================================================
preflight_and_submit() {
  autodetect
  load_tools
  [[ -n "$REGENIE_CMD" ]] || die "regenie not found — run ./run_gwas.sh --check for details"

  #--- required inputs -------------------------------------------------
  # Every path is made absolute: SLURM stages run from a different cwd,
  # so a relative --pheno/--covar/--keep would break inside the job.
  [[ -n "$PHENO" ]] || { usage; echo; die "--pheno is required"; }
  (( KEEP_SET )) || KEEP="$KEEP_DEFAULT"

  need_file "$PHENO" "Phenotype file"
  PHENO="$(abspath "$PHENO")"
  if [[ -n "$COVAR" ]]; then need_file "$COVAR" "Covariate file"; COVAR="$(abspath "$COVAR")"; fi
  if [[ -n "$KEEP"  ]]; then need_file "$KEEP"  "Sample keep list"; KEEP="$(abspath "$KEEP")"; fi

  set_geno_args "$GENO_STEP1" "step 1"     # dies with a clear message if absent
  set_geno_args "$GENO_STEP2" "step 2"

  if [[ ! -f "$PRUNED_SNPS" && -z "$PLINK2_CMD" ]]; then
    die "no cached LD-pruned SNP list at ${PRUNED_SNPS}, and PLINK2 was not found.
     Pass --plink2 /path/to/plink2 to build one (a one-off cost, cached for
     every later phenotype), or --pruned-snps to reuse a labmate's list."
  fi

  #--- phenotype header ------------------------------------------------
  local ncol traits_avail
  ncol=$(awk 'NR==1 { print NF; exit }' "$PHENO")
  traits_avail=$(awk 'NR==1 { for (i=3; i<=NF; i++) printf "%s ", $i; exit }' "$PHENO")
  awk 'NR==1 { if ($1 != "FID" || $2 != "IID") exit 1; exit 0 }' "$PHENO" \
    || die "phenotype file must start with columns FID IID — found: $(head -1 "$PHENO")"

  #--- resolve the trait column ----------------------------------------
  if [[ -z "$TRAIT" ]]; then
    (( ncol == 3 )) || die "--trait is required: ${PHENO} has ${ncol} columns.
     Available: ${traits_avail}"
    TRAIT=$(awk 'NR==1 { print $3; exit }' "$PHENO")
    say "auto-detected phenotype column: ${TRAIT}"
  fi
  [[ "$TRAIT" =~ ^[A-Za-z0-9_.-]+$ ]] \
    || die "trait name '${TRAIT}' would break filenames — use A-Z a-z 0-9 _ . - only"

  local tcol
  tcol=$(awk -v t="$TRAIT" 'NR==1 { for (i=1; i<=NF; i++) if ($i == t) { print i; exit } }' "$PHENO")
  [[ -n "$tcol" ]] || die "column '${TRAIT}' not found in ${PHENO}
     Available: ${traits_avail}"

  #--- distribution + sample overlap, in one pass ----------------------
  local stats
  stats=$(awk -v c="$tcol" -v keepfile="${KEEP}" '
    BEGIN {
      usekeep = (keepfile != "")
      if (usekeep) while ((getline line < keepfile) > 0) { split(line, a, /[ \t]+/); K[a[2]] = 1 }
    }
    NR > 1 {
      total++
      v = $c
      if (v == "NA" || v == "" || v == "-9")            { miss++; next }
      if (v !~ /^-?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) { bad++;  next }
      n++; s += v; ss += v * v
      if (n == 1 || v + 0 < min) min = v + 0
      if (n == 1 || v + 0 > max) max = v + 0
      if (v + 0 != 0 && v + 0 != 1) notbin01 = 1
      if (v + 0 != 1 && v + 0 != 2) notbin12 = 1
      if (!(v in seen) && nu < 20) { seen[v] = 1; nu++ }
      if (!usekeep || ($2 in K)) inkeep++
    }
    END {
      mean = (n ? s / n : 0)
      sd   = (n > 1 ? sqrt((ss - n * mean * mean) / (n - 1)) : 0)
      printf "%d %d %d %d %d %.6g %.6g %.6g %.6g %d %d %d",
             total, n, miss, bad, inkeep, mean, sd, min, max,
             (notbin01 ? 0 : 1), (notbin12 ? 0 : 1), nu
    }' "$PHENO")
  local P_TOTAL P_N P_MISS P_BAD P_INKEEP P_MEAN P_SD P_MIN P_MAX IS_BIN01 IS_BIN12 P_NUNIQ
  read -r P_TOTAL P_N P_MISS P_BAD P_INKEEP P_MEAN P_SD P_MIN P_MAX IS_BIN01 IS_BIN12 P_NUNIQ <<<"$stats"

  (( P_N > 0 ))    || die "no numeric values found in column '${TRAIT}'"
  (( P_BAD == 0 )) || warn "${P_BAD} non-numeric values in '${TRAIT}' will be dropped by REGENIE"

  if (( IS_BIN01 == 1 || IS_BIN12 == 1 )) && (( FORCE_BINARY == 0 )); then
    die "'${TRAIT}' looks BINARY (only $( ((IS_BIN01)) && echo '0/1' || echo '1/2' ) values present).
     This script runs a linear-model quantitative GWAS. A binary trait needs
     REGENIE --bt with Firth correction. Use --force-binary to override."
  fi

  if [[ -n "$KEEP" ]]; then
    (( P_INKEEP > 0 )) || die "no phenotyped sample appears in the keep list ${KEEP}
     Check that column 2 (IID) of both files holds the same UKB EIDs."
  fi
  (( P_INKEEP >= 1000 )) || warn "only ${P_INKEEP} phenotyped samples — badly underpowered for GWAS"

  #--- covariate sanity ------------------------------------------------
  if [[ -n "$COVAR" ]]; then
    awk 'NR==1 { if ($1 != "FID" || $2 != "IID") exit 1; exit 0 }' "$COVAR" \
      || die "covariate file must start with columns FID IID: ${COVAR}"
  else
    warn "no covariate file — nothing adjusts for ancestry PCs, age or sex.
     Population structure will inflate your results. Pass --covar unless the
     phenotype is already fully residualised."
  fi

  #--- output directory ------------------------------------------------
  OUTDIR="${OUTDIR:-${OUT_ROOT}/${TRAIT}}"
  mkdir -p "${OUTDIR}/logs"
  OUTDIR="$(cd "$OUTDIR" && pwd)"

  #--- report ----------------------------------------------------------
  cat <<EOF

  ---------------------------------------------------------------
   GWAS: ${TRAIT}
  ---------------------------------------------------------------
   phenotype file    ${PHENO}
   column            ${TRAIT}  (column ${tcol} of ${ncol})
   non-missing       ${P_N} of ${P_TOTAL} rows   (${P_MISS} NA)
   $( [[ -n "$KEEP" ]] && echo "overlap with keep ${P_INKEEP} samples  <- expected GWAS N" \
                       || echo "expected GWAS N   ${P_INKEEP} (limited by genotyped samples)" )
   distribution      mean ${P_MEAN}  sd ${P_SD}  range [${P_MIN}, ${P_MAX}]
   distinct values   $( (( P_NUNIQ >= 20 )) && echo "20+" || echo "${P_NUNIQ}" )

   covariates        ${COVAR:-<none>}
                     ${COVAR_COLS:-<all columns>}$( [[ -n "$CAT_COVAR" ]] && echo "   categorical: ${CAT_COVAR}" )
   keep list         ${KEEP:-<all genotyped samples>}
   genotypes step 2  $(basename "${GENO_STEP2}")  [$(geno_format "$GENO_STEP2")]
   genotypes step 1  $(basename "${GENO_STEP1}")  [$(geno_format "$GENO_STEP1")]
   pruned SNP list   $( [[ -f "$PRUNED_SNPS" ]] && echo "${PRUNED_SNPS} ($(wc -l < "$PRUNED_SNPS") SNPs)" \
                                                || echo "${PRUNED_SNPS} (will be built in step 1)" )
   chromosomes       ${CHR_FIRST}-${CHR_LAST}
   rank-inv-normal   $( ((APPLY_RINT)) && echo yes || echo "no (--no-rint)" )
   min MAC           ${MIN_MAC}$( [[ -n "$MAF" ]] && echo "   MAF floor ${MAF}" )

   output            ${OUTDIR}
   scheduler         ${ACCOUNT:+account ${ACCOUNT}  }partition ${PARTITION}
  ---------------------------------------------------------------

EOF

  #--- freeze the configuration ----------------------------------------
  local cfg="${OUTDIR}/run_config.sh"
  cat > "$cfg" <<EOF
# run_gwas.sh ${GWAS_VERSION} — resolved configuration for ${TRAIT}
# written $(date '+%Y-%m-%d %H:%M:%S') by ${USER:-unknown}
# rerun a stage with:
#   ${SCRIPT_PATH} --stage <step1|step2|merge> --config ${cfg}
SCRIPT_PATH="${SCRIPT_PATH}"
PHENO="${PHENO}"
TRAIT="${TRAIT}"
OUTDIR="${OUTDIR}"
COVAR="${COVAR}"
COVAR_COLS="${COVAR_COLS}"
CAT_COVAR="${CAT_COVAR}"
KEEP="${KEEP}"
GENO_STEP1="${GENO_STEP1}"
GENO_STEP2="${GENO_STEP2}"
SAMPLE_FILE="${SAMPLE_FILE}"
PRUNED_SNPS="${PRUNED_SNPS}"
CHROMS="${CHROMS}"
MAF="${MAF}"
MIN_MAC="${MIN_MAC}"
APPLY_RINT=${APPLY_RINT}
GW_SIG="${GW_SIG}"
CLUMP_KB=${CLUMP_KB}
BSIZE_STEP1=${BSIZE_STEP1}
BSIZE_STEP2=${BSIZE_STEP2}
ACCOUNT="${ACCOUNT}"
PARTITION="${PARTITION}"
TIME_STEP2="${TIME_STEP2}"
NTASK_STEP2=${NTASK_STEP2}
THREADS="${THREADS}"
CONDA_ENV="${CONDA_ENV}"
REGENIE_BIN="${REGENIE_BIN}"
PLINK2_BIN="${PLINK2_BIN}"
SCRATCH_ROOT="${SCRATCH_ROOT}"
OUT_ROOT="${OUT_ROOT}"
EOF

  #--- build the submission chain --------------------------------------
  # sbatch options must precede the script path, so they live in their
  # own arrays and are composed at call time.
  local L="${OUTDIR}/logs"
  local -a acct=() mail_fail=() mail_done=()
  [[ -n "$ACCOUNT" ]] && acct=(-A "$ACCOUNT")
  if [[ -n "$EMAIL" ]]; then
    mail_fail=(--mail-type=FAIL --mail-user="$EMAIL")
    mail_done=(--mail-type=END,FAIL --mail-user="$EMAIL")
  fi

  local -a opt1=("${acct[@]+"${acct[@]}"}" -p "$PARTITION" -N 1 -n "$NTASK_STEP1"
                 -t "$TIME_STEP1" -J "${TRAIT}_s1"
                 -o "${L}/step1_%j.out" -e "${L}/step1_%j.err"
                 "${mail_fail[@]+"${mail_fail[@]}"}")
  local -a opt2=(--array="${CHR_FIRST}-${CHR_LAST}" "${acct[@]+"${acct[@]}"}" -p "$PARTITION"
                 -N 1 -n "$NTASK_STEP2" -t "$TIME_STEP2" -J "${TRAIT}_s2"
                 -o "${L}/step2_%A_%a.out" -e "${L}/step2_%A_%a.err"
                 "${mail_fail[@]+"${mail_fail[@]}"}")
  local -a opt3=("${acct[@]+"${acct[@]}"}" -p "$PARTITION" -N 1 -n "$NTASK_MERGE"
                 -t "$TIME_MERGE" -J "${TRAIT}_mg"
                 -o "${L}/merge_%j.out" -e "${L}/merge_%j.err"
                 "${mail_done[@]+"${mail_done[@]}"}")

  local -a arg1=("$SCRIPT_PATH" --stage step1 --config "$cfg")
  local -a arg2=("$SCRIPT_PATH" --stage step2 --config "$cfg")
  local -a arg3=("$SCRIPT_PATH" --stage merge --config "$cfg")

  local step1_done=0
  [[ -f "${OUTDIR}/step1/step1_${TRAIT}_pred.list" ]] && step1_done=1

  if (( DRY_RUN )); then
    say "DRY RUN — inputs validated, nothing submitted."
    echo
    if (( step1_done )); then
      echo "  step 1 : already complete, would be skipped"
    else
      echo "  step 1 : sbatch ${opt1[*]} ${arg1[*]}"
    fi
    echo "  step 2 : sbatch ${opt2[*]} --dependency=afterok:<step1> ${arg2[*]}"
    echo "  merge  : sbatch ${opt3[*]} --dependency=afterok:<step2> ${arg3[*]}"
    echo
    return 0
  fi

  local jid1="" jid2="" jid3=""
  if (( step1_done )); then
    say "step 1 predictions already present — skipping to step 2"
    jid2=$(submit_job "${opt2[@]}" "${arg2[@]}")
  else
    jid1=$(submit_job "${opt1[@]}" "${arg1[@]}")
    say "submitted step 1   job ${jid1}"
    jid2=$(submit_job "${opt2[@]}" --dependency="afterok:${jid1}" "${arg2[@]}")
  fi
  say "submitted step 2   job ${jid2}  (array ${CHR_FIRST}-${CHR_LAST})"

  jid3=$(submit_job "${opt3[@]}" --dependency="afterok:${jid2}" "${arg3[@]}")
  say "submitted merge    job ${jid3}"

  cat <<EOF

  Chain submitted — nothing else to do. Results will appear in:
      ${OUTDIR}/${TRAIT}.gwas.tsv.gz
      ${OUTDIR}/${TRAIT}.tophits.tsv
      ${OUTDIR}/${TRAIT}.summary.txt

  Watch progress:   squeue -u \$USER
  Cancel this run:  scancel ${jid1:+$jid1 }${jid2} ${jid3}
  Logs:             ${OUTDIR}/logs/

EOF
}

#=======================================================================
#  DISPATCH
#=======================================================================
[[ "$MODE" == "check" ]] && { run_check; exit $?; }

case "$STAGE" in
  all)   preflight_and_submit ;;
  step1) [[ -n "$CONFIG" ]] || die "--stage step1 requires --config <run_config.sh>"; run_step1 ;;
  step2) [[ -n "$CONFIG" ]] || die "--stage step2 requires --config <run_config.sh>"; run_step2 ;;
  merge) [[ -n "$CONFIG" ]] || die "--stage merge requires --config <run_config.sh>"; run_merge ;;
  *)     die "unknown --stage '${STAGE}' — use all | step1 | step2 | merge" ;;
esac
