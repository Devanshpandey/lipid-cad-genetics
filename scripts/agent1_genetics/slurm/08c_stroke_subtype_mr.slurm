#!/bin/bash
#SBATCH -J strokemr
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/strokemr_%j.out
#SBATCH -e /path/to/cad-genetics/strokemr_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
M=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
D=$_REPO/data/external/megastroke
mkdir -p $D $_REPO/exports/agent1/stroke

# MEGASTROKE European-ancestry summary statistics (Malik et al., Nat Genet 2018;
# PMID 29531354), downloaded from the GWAS Catalog FTP. The accession-to-subtype
# mapping is NOT sequential -- it is as listed below.
FTP=ftp://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000
declare -A ACC=( [MEGASTROKE.1.AS.EUR.out]=GCST006906    # any stroke
                 [MEGASTROKE.2.AIS.EUR.out]=GCST006908   # any ischemic stroke
                 [MEGASTROKE.3.LAS.EUR.out]=GCST006907   # large-artery atherosclerotic
                 [MEGASTROKE.4.CES.EUR.out]=GCST006910   # cardioembolic
                 [MEGASTROKE.5.SVS.EUR.out]=GCST006909 ) # small-vessel
for F in "${!ACC[@]}"; do
  [[ -s $D/$F ]] || curl -sS --retry 3 -o $D/$F "$FTP/${ACC[$F]}/$F"
done
# LDL-C is the positive control: its association with large-artery atherosclerotic
# stroke is established, so if the subtype data are informative LDL-C should show it.
declare -A S=( [AS]=MEGASTROKE.1.AS.EUR.out [AIS]=MEGASTROKE.2.AIS.EUR.out \
               [LAS]=MEGASTROKE.3.LAS.EUR.out [CES]=MEGASTROKE.4.CES.EUR.out \
               [SVS]=MEGASTROKE.5.SVS.EUR.out )
for K in AS AIS LAS CES SVS; do
  echo "=== $K ==="
  Rscript $_REPO/scripts/agent1_genetics/R/run_external_mr.R \
    --sumstat_dir $M --outcome_file $D/${S[$K]} \
    --outcome_format megastroke --outcome_label MEGASTROKE_$K \
    --exposures LDL_C,LPA --mode univariable \
    --ld_ref $SCRATCH/lpa_ldref --plink /path/to/software/plink \
    --cis_dense "LPA:6:160500000:161700000:0.01:1000" \
    --out $_REPO/exports/agent1/stroke/mr_${K}.csv
done
echo ALL_DONE
