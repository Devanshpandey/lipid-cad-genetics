#!/bin/bash
set -euo pipefail
ROOT=/path/to/cad-genetics
PH=$ROOT/results/agent1_genetics/phenotypes
OUT=$ROOT/results/agent1_genetics/strengthen
PLINK2=/path/to/software/plink2
RSCRIPT=/path/to/anaconda3/envs/r_env/bin/Rscript
WES=/corral-repl/utexas/UKB-Imaging-Genetics/UKB_GENOTYPE/ukb_pfile/wes_allchr
mkdir -p $OUT; cd $OUT
awk '$5=="PDE3B" && $6=="pLoF"{print $1":"$2":"$3":"$4}' \
  $ROOT/results/agent1_genetics/burden_masks/by_chr/anno_chr11.txt | sort -u > pde3b_plof.snplist
echo "pLoF variants to extract: $(wc -l < pde3b_plof.snplist)"
$PLINK2 --pfile $WES --extract pde3b_plof.snplist --export A --out pde3b_plof 2>&1 | tail -4
echo "variants present in .raw: $(( $(head -1 pde3b_plof.raw | tr '\t' '\n' | wc -l) - 6 ))"
$RSCRIPT $OUT/pde3b_phewas.R
