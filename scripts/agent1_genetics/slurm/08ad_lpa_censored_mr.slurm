#!/bin/bash
#------------------------------------------------------------
# Rerun the published Lp(a) analyses on the censored-block GWAS from 08ac, with the published scripts unchanged.
# Arms, each a --sumstat_dir differing from merged_statinfix in lipids_LPA only:
#   ctrl     merged_lpameasctrl       LPA_MEAS; must reproduce the published numbers
#   fixediv  merged_lpacens_fixediv   published instruments, corrected BETA/SE: the pure effect of the correction
#   cens     merged_lpacens           instruments reselected on the censored GWAS: what would ship
# Commands mirror 08h (MVMR, FinnGen, CARDIoGRAM, MEGASTROKE, covariance), 08p (subtype contrasts), 08m (cis-MR)
# and 08e (scores). Outputs go to exports/agent1/lpa_censored/<arm>/; published exports are not touched.
#   sbatch --dependency=afterok:<08ac stage B> 08ad_lpa_censored_mr.slurm
#------------------------------------------------------------
#SBATCH -J lpacens_mr
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 128
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/logs/lpacens_mr_%j.out
#SBATCH -e /path/to/cad-genetics/logs/lpacens_mr_%j.err
set -uo pipefail
umask 077
_REPO=/path/to/cad-genetics
source "${_REPO}/config/tacc_paths.sh"
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
export R_DATATABLE_NUM_THREADS=8 OMP_NUM_THREADS=8
R=$_REPO/scripts/agent1_genetics/R
G=$_REPO/results/agent1_genetics/gwas_summary_stats
PH=$_REPO/results/agent1_genetics/phenotypes
E=$_REPO/exports/agent1/lpa_censored
PLINK=/path/to/software/plink
LDREF=$SCRATCH/lpa_ldref
CISREF=$SCRATCH/cis_ldref
ST=$SCRATCH/statint
RAW=/path/to/ukb_phenotypes
FG=$_REPO/data/external/finngen/finngen_R11_I9_CHD.gz
C4D=$_REPO/data/external/cardiogram/cad.add.160614.website.txt
MS=$_REPO/data/external/megastroke
CIS='LPA:6:160500000:161700000:0.01:1000'
JOBS=$E/jobs.txt; : > "$JOBS"

declare -A DIR=( [ctrl]=$G/merged_lpameasctrl [fixediv]=$G/merged_lpacens_fixediv [cens]=$G/merged_lpacens )
declare -A MSF=( [AS]=MEGASTROKE.1.AS.EUR.out [AIS]=MEGASTROKE.2.AIS.EUR.out [LAS]=MEGASTROKE.3.LAS.EUR.out
                 [CES]=MEGASTROKE.4.CES.EUR.out [SVS]=MEGASTROKE.5.SVS.EUR.out )
add() { local log="$1"; shift; echo "$* > $E/logs/$log.log 2>&1 || echo FAILED $log" >> "$JOBS"; }

for ARM in ctrl fixediv cens; do
  S=${DIR[$ARM]}; O=$E/$ARM; mkdir -p "$O/stroke" "$E/logs"
  [ -s "$S/lipids_LPA.txt.gz" ] || { echo "MISSING $S/lipids_LPA.txt.gz"; exit 1; }
  if [ "$ARM" = ctrl ]; then OUTS="CAD"; else OUTS="CAD MI REVASC STROKE HF CV_DEATH MACE PAD"; fi
  for T in $OUTS; do
    add mvmr_${ARM}_$T Rscript $R/run_mvmr.R --outcome $T --sumstat_dir $S --out $O/mvmr_$T.csv
  done
  add finngen_$ARM Rscript $R/run_external_mr.R --sumstat_dir $S --outcome_file $FG --outcome_format finngen \
    --outcome_label FinnGen_R11_CAD --exposures LDL_C,HDL_C,TRIGLY,LPA --mode mvmr --ld_ref $LDREF --plink $PLINK \
    --cis_dense "$CIS" --out $O/mvmr_finngen.csv
  add c4d_$ARM Rscript $R/run_external_mr.R --sumstat_dir $S --outcome_file $C4D --outcome_format cardiogram \
    --outcome_label CARDIoGRAMplusC4D --exposures LPA --mode univariable --ld_ref $LDREF --plink $PLINK \
    --cis_dense "$CIS" --out $O/mr_cardiogram.csv
  add contrasts_$ARM Rscript $R/stroke_subtype_contrasts.R --sumstat_dir $S --megastroke_dir $MS --ld_ref $CISREF \
    --plink $PLINK --nboot 1000 --out_prefix $O/stroke_subtype_contrasts
  if [ "$ARM" = ctrl ]; then P=$PH/pheno_quantitative_statinfix.txt; else P=$PH/lpa_censored/pheno_quantitative_statinfix_lpacens.txt; fi
  add cov_$ARM Rscript $R/mvmr_covariance.R --sumstat_dir $S --pheno $P --outcomes CAD,STROKE,MI,REVASC --nboot 1000 \
    --cis_list $ST/lpa_cis_ids.txt --out $O/mvmr_covariance.csv
  if [ "$ARM" != ctrl ]; then
    for K in AS AIS LAS CES SVS; do
      add stroke_${ARM}_$K Rscript $R/run_external_mr.R --sumstat_dir $S --outcome_file $MS/${MSF[$K]} --outcome_format megastroke \
        --outcome_label MEGASTROKE_$K --exposures LPA --mode univariable --ld_ref $LDREF --plink $PLINK \
        --cis_dense "$CIS" --out $O/stroke/mr_$K.csv
    done
    add cismr_$ARM Rscript $R/run_cis_mr.R --sumstat_dir $S --outcome_dir $S --ldl_file lipids_LDL_C.txt.gz \
      --ld_ref $LDREF --plink $PLINK --label $ARM --out $O/cis_mr.csv
  fi
done

# Statin score: the published 45 cis variants with censored-GWAS weights. ctrl keeps the published weights file.
Rscript -e 'suppressPackageStartupMessages(library(data.table)); a <- commandArgs(TRUE)
  w <- fread(a[1], header = FALSE, col.names = c("SNP", "EA", "B")); s <- fread(cmd = paste("gunzip -c", shQuote(a[2])), select = c("ID", "ALLELE1", "BETA"))
  o <- s[ID %in% w$SNP, .(ID, ALLELE1, BETA)]; stopifnot(nrow(o) == nrow(w), !anyNA(o$BETA)); fwrite(o, a[3], sep = "\t", col.names = FALSE)' \
  "$ST/lpa_score.txt" "${DIR[cens]}/lipids_LPA.txt.gz" "$E/cens/lpa_score_cens.txt"
for ARM in ctrl cens; do
  if [ "$ARM" = ctrl ]; then W=$ST/lpa_score.txt; else W=$E/cens/lpa_score_cens.txt; fi
  add statin_$ARM Rscript $R/statin_interaction_incident.R --scores $ST/dose_LPA.raw,$ST/dose_LDL.raw --weights $W,$ST/ldl_score.txt \
    --raw_dir $RAW --covar $PH/covariates.txt --out $E/$ARM/statin_incident.csv
done

echo "[$(date '+%H:%M:%S')] $(wc -l < "$JOBS") jobs"
xargs -d '\n' -P 12 -I{} bash -c '{}' < "$JOBS"
echo "[$(date '+%H:%M:%S')] finished; failures, if any, are printed above"
echo ALL_DONE
