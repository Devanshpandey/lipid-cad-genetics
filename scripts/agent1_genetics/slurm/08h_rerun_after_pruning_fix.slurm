#!/bin/bash
#------------------------------------------------------------
# Re-run every Mendelian randomization analysis after the distance-pruning fix.
#
# The guard `any(same_chr & abs(keep_pos[same_chr] - pos) < window)` compared two
# vectors of different lengths, so R recycled the shorter one and the guard could
# return FALSE for variants that were in fact inside the window. Correlated
# variants therefore entered instrument sets that are assumed independent.
# Every distance-pruned analysis in the project is affected.
#
# Results go to *_prunefix paths so the pre-fix outputs stay on disk for comparison.
#------------------------------------------------------------
#SBATCH -J prunefix
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 08:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/prunefix_%j.out
#SBATCH -e /path/to/cad-genetics/prunefix_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
E=$_REPO/exports/agent1
PLINK=/path/to/software/plink
LDREF=$SCRATCH/lpa_ldref
mkdir -p $E/mvmr_prunefix $E/mvmr_apob_prunefix $E/stroke_prunefix

echo "############ [1] primary four-exposure MVMR, 8 endpoints"
for T in CAD MI REVASC STROKE HF CV_DEATH MACE PAD; do
  echo "=== $T ==="
  Rscript $_REPO/scripts/agent1_genetics/R/run_mvmr.R --outcome $T --sumstat_dir $S \
    --out $E/mvmr_prunefix/mvmr_$T.csv
done

echo "############ [2] focused LDL-C / ApoB / triglyceride model"
ln -sf $_REPO/results/agent1_genetics/gwas_summary_stats/merged/lipids_APOB.txt.gz $S/lipids_APOB.txt.gz
for T in CAD MI MACE; do
  echo "=== $T ==="
  Rscript $_REPO/scripts/agent1_genetics/R/run_mvmr_apob.R --outcome $T --sumstat_dir $S \
    --out $E/mvmr_apob_prunefix/mvmr_apob_$T.csv
done

echo "############ [3] FinnGen R11 two-sample MVMR"
Rscript $_REPO/scripts/agent1_genetics/R/run_external_mr.R \
  --sumstat_dir $S --outcome_file $_REPO/data/external/finngen/finngen_R11_I9_CHD.gz \
  --outcome_format finngen --outcome_label FinnGen_R11_CAD \
  --exposures LDL_C,HDL_C,TRIGLY,LPA --mode mvmr \
  --ld_ref $LDREF --plink $PLINK \
  --cis_dense "LPA:6:160500000:161700000:0.01:1000" \
  --out $E/mvmr_finngen_prunefix.csv

echo "############ [4] CARDIoGRAMplusC4D univariable"
Rscript $_REPO/scripts/agent1_genetics/R/run_external_mr.R \
  --sumstat_dir $_REPO/results/agent1_genetics/gwas_summary_stats/merged \
  --outcome_file $_REPO/data/external/cardiogram/cad.add.160614.website.txt \
  --outcome_format cardiogram --outcome_label CARDIoGRAMplusC4D \
  --exposures LDL_C,HDL_C,TRIGLY,LPA --mode univariable \
  --ld_ref $LDREF --plink $PLINK \
  --cis_dense "LPA:6:160500000:161700000:0.01:1000" \
  --out $E/mr_cardiogram_prunefix.csv

echo "############ [5] MEGASTROKE subtypes"
D=$_REPO/data/external/megastroke
declare -A MS=( [AS]=MEGASTROKE.1.AS.EUR.out [AIS]=MEGASTROKE.2.AIS.EUR.out \
                [LAS]=MEGASTROKE.3.LAS.EUR.out [CES]=MEGASTROKE.4.CES.EUR.out \
                [SVS]=MEGASTROKE.5.SVS.EUR.out )
for K in AS AIS LAS CES SVS; do
  echo "=== $K ==="
  Rscript $_REPO/scripts/agent1_genetics/R/run_external_mr.R \
    --sumstat_dir $S --outcome_file $D/${MS[$K]} \
    --outcome_format megastroke --outcome_label MEGASTROKE_$K \
    --exposures LDL_C,LPA --mode univariable \
    --ld_ref $LDREF --plink $PLINK \
    --cis_dense "LPA:6:160500000:161700000:0.01:1000" \
    --out $E/stroke_prunefix/mr_$K.csv
done

echo "############ [6] covariance sensitivity"
Rscript $_REPO/scripts/agent1_genetics/R/mvmr_covariance.R \
  --sumstat_dir $S \
  --pheno $_REPO/results/agent1_genetics/phenotypes/pheno_quantitative_statinfix.txt \
  --outcomes CAD,STROKE,MI,REVASC \
  --nboot 1000 \
  --cis_list /path/to/scratch/statint/lpa_cis_ids.txt \
  --out $E/mvmr_covariance_prunefix.csv

echo "############ [7] LD audit of the corrected instrument set"
awk 'NR>1{print $1}' $E/mvmr_prunefix/mvmr_CAD.instruments.tsv > $SCRATCH/ivs_fixed.txt
$PLINK --bfile $LDREF --extract $SCRATCH/ivs_fixed.txt --r2 inter-chr --ld-window-r2 0.005 \
       --out $SCRATCH/ivld_fixed --threads 16 --memory 40000 > /dev/null 2>&1 || true
echo "instruments: $(wc -l < $SCRATCH/ivs_fixed.txt)"
if [ -f $SCRATCH/ivld_fixed.ld ]; then
  echo "pairs with r2>0.005: $(( $(wc -l < $SCRATCH/ivld_fixed.ld) - 1 ))"
  sort -k7,7gr $SCRATCH/ivld_fixed.ld | head -5
else
  echo "no pair exceeded r2 0.005"
fi
echo ALL_DONE
