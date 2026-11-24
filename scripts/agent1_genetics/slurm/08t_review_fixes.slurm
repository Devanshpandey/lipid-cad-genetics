#!/bin/bash
#SBATCH -J revfix
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 00:45:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/revfix_%j.out
#SBATCH -e /path/to/cad-genetics/revfix_%j.err
# ============================================================================
# Two numbers peer review showed the manuscript states inconsistently.
#
# (1) TABLE S1 STATIN COUNT. Table S1 reports 27,549 statin users, which is the
#     superseded field-6153-only count (UK Biobank asks 6153 of women only).
#     Every analysis now uses 6153 U 6177. Recount over the ANALYSED cohort.
#
# (2) MAX PAIRWISE r2 OF THE cis-LPA INSTRUMENT. Methods say 0.0098, Table S39
#     says 0.0115, for what is described as the same instrument. Recompute for
#     each set explicitly so the manuscript can state the right one -- or state
#     both, if they are genuinely different sets.
# ============================================================================
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
PLINK=/path/to/software/plink
S=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
CV=$_REPO/results/agent1_genetics/phenotypes/covariates.txt
W=$SCRATCH/revfix; mkdir -p $W

echo "########## (1) STATIN COUNT, corrected union (6153 U 6177) ##########"
awk 'NR>1{print $2}' $CV | sort -u > $W/analysed.ids
echo "analysed cohort (covariates.txt): $(wc -l < $W/analysed.ids)"
awk 'NR==FNR{keep[$1];next} FNR>1 && ($1 in keep){n++; if($2==1)s++} END{
  printf "  in-cohort with statin field: %d\n  STATIN USERS (6153 U 6177): %d (%.1f%%)\n", n, s, 100*s/n}' \
  $W/analysed.ids /path/to/scratch/statint/statin.tsv

echo
echo "########## (2) MAX PAIRWISE r2 OF EACH cis-LPA INSTRUMENT SET ##########"
REF=/path/to/scratch/lpa_ldref
gunzip -c $S/lipids_LPA.txt.gz | awk 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
  $h["CHROM"]==6 && $h["GENPOS"]>=160500000 && $h["GENPOS"]<=161700000 &&
  $h["LOG10P"]>7.30103 && !($h["BETA"]=="NA"){print $h["ID"], $h["LOG10P"]}' > $W/c.raw
awk 'NR==FNR{if(FNR==1){mx=$2;mn=$2} if($2>mx)mx=$2; if($2<mn)mn=$2; next}
  {p=10^(-(1+($2-mn)/(mx-mn)*299)); print $1,p}' $W/c.raw $W/c.raw > $W/c.body
(echo "SNP P"; cat $W/c.body) > $W/c.assoc

$PLINK --bfile $REF --clump $W/c.assoc --clump-p1 1 --clump-p2 1 \
       --clump-r2 0.01 --clump-kb 1000 --out $W/cl01 >/dev/null 2>&1
awk 'NR>1&&NF{print $3}' $W/cl01.clumped > $W/set45.txt
echo "clumped at r2<0.01: $(wc -l < $W/set45.txt) variants"

maxr2 () {  # $1 = id list, $2 = label
  local n=$(wc -l < "$1")
  if [ "$n" -lt 2 ]; then echo "  $2: n=$n, no pairs"; return; fi
  $PLINK --bfile $REF --extract "$1" --r2 --ld-window 99999 --ld-window-kb 2000 \
         --ld-window-r2 0 --out "$W/r2_$2" >/dev/null 2>&1
  awk -v L="$2" -v N="$n" 'NR>1&&NF{if($7>m)m=$7; c++} END{
    printf "  %-14s n=%-4s pairs=%-7s max_r2=%.5f\n", L, N, c+0, m}' "$W/r2_$2.ld"
}
maxr2 $W/set45.txt set45

# the 23-variant contrast set: present AND non-missing in all five subtype files
MS=$_REPO/data/external/megastroke
cp $W/set45.txt $W/common.txt
for P in LAS CES SVS AIS AS; do
  f=$(ls $MS/*${P}* 2>/dev/null | head -1)
  [ -z "$f" ] && { echo "  (no MEGASTROKE file for $P, skipping)"; continue; }
  if [[ "$f" == *.gz ]]; then RD="gunzip -c $f"; else RD="cat $f"; fi
  $RD | awk 'NR>1 && $0!~/NA/{print $1}' | sort -u > $W/ids_$P
  comm -12 <(sort -u $W/common.txt) $W/ids_$P > $W/tmp && mv $W/tmp $W/common.txt
done
echo "common across all five subtype files: $(wc -l < $W/common.txt) variants"
maxr2 $W/common.txt common
echo ALL_DONE
