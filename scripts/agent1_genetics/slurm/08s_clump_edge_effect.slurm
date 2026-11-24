#!/bin/bash
#SBATCH -J clumpedge
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 01:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/clumpedge_%j.out
#SBATCH -e /path/to/cad-genetics/clumpedge_%j.err
# ============================================================================
# Does the non-monotone clumping count come from the clump RADIUS being smaller
# than the WINDOW being clumped?
#
# Observed: in the 1.2 Mb subtype-contrast window, r2<0.01 -> 99 instruments but
# r2<0.001 -> 117. Greedy clumping at a lower r2 absorbs a strict superset per
# clump, so it should leave FEWER index variants, never more.
#
# HYPOTHESIS. --clump-kb 1000 means a clump reaches only 1 Mb from its index. In a
# 1.2 Mb window a variant can therefore sit beyond every index's reach and survive
# regardless of LD. Which variants become indices changes with the r2 threshold, so
# the set of out-of-reach variants changes too, and the count can move either way.
#
# PREDICTION. Widen the radius past the window span (>=1200 kb) and monotonicity is
# restored. The 335 kb cis-MR window already sits inside a single 1 Mb radius, so it
# should be monotone at EVERY radius -- that is the control arm. If the contrast stays
# non-monotone at 5000 kb, or the cis-MR window is non-monotone anywhere, the
# hypothesis is wrong and the cis-MR counts need re-examining too.
# ============================================================================
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
PLINK=/path/to/software/plink
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
W=$SCRATCH/clumpedge; mkdir -p $W
REF=$SCRATCH/cis_ldref_full

gunzip -c $S/lipids_LPA.txt.gz | awk 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $h["CHROM"]==6 && $h["GENPOS"]>=160852515 && $h["GENPOS"]<=161187407 && $h["LOG10P"]>7.30103 && $h["ID"]~/^rs/ {p=10^(-$h["LOG10P"]); if(p==0)p=1e-300; print $h["ID"], p}' > $W/b1
(echo "SNP P"; cat $W/b1) > $W/cismr.assoc
gunzip -c $S/lipids_LPA.txt.gz | awk 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $h["CHROM"]==6 && $h["GENPOS"]>=160500000 && $h["GENPOS"]<=161700000 && $h["LOG10P"]>7.30103 && !($h["BETA"]=="NA") {print $h["ID"], $h["LOG10P"]}' > $W/c.raw
awk 'NR==FNR{if(FNR==1){mx=$2;mn=$2} if($2>mx)mx=$2; if($2<mn)mn=$2; next}{p=10^(-(1+($2-mn)/(mx-mn)*299)); print $1, p}' $W/c.raw $W/c.raw > $W/c.body
(echo "SNP P"; cat $W/c.body) > $W/contrast.assoc
echo "candidates: cis-MR(335kb)=$(wc -l < $W/b1)  contrast(1.2Mb)=$(wc -l < $W/c.body)"

cnt(){ $PLINK --bfile $REF --clump $1 --clump-p1 1 --clump-r2 $2 --clump-kb $3 --out $W/sw >/dev/null 2>&1; awk 'NR>1&&NF{n++}END{print n+0}' $W/sw.clumped; }
for SET in contrast cismr; do
  [ $SET = contrast ] && SPAN="1.2 Mb span" || SPAN="335 kb span"
  echo "=== $SET window ($SPAN), clump radius swept ==="
  for KB in 400 600 1000 1200 2000 5000; do
    a=$(cnt $W/$SET.assoc 0.01 $KB); b=$(cnt $W/$SET.assoc 0.001 $KB)
    if [ "$b" -le "$a" ]; then V="monotone OK"; else V="NON-MONOTONE"; fi
    printf "  radius %-5s kb :  r2<0.01 = %-4s  r2<0.001 = %-4s  %s\n" "$KB" "$a" "$b" "$V"
  done
done
echo ALL_DONE
