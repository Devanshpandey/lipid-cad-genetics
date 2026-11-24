#!/bin/bash
#SBATCH -J statint
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 32
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/statint_%j.out
#SBATCH -e /path/to/cad-genetics/statint_%j.err
set -euo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
M=$_REPO/results/agent1_genetics/gwas_summary_stats/merged_statinfix
W=$SCRATCH/statint; mkdir -p $W; export W

echo "[1] instrument + weight files"
# Lp(a): locus-dense cis instrument, same construction as the MVMR (r2<0.01 within
# 1 Mb of LPA). Ranked on LOG10P, never on pval -- pval underflows here.
zcat $M/lipids_LPA.txt.gz | awk 'NR>1 && $1==6 && $2>=160500000 && $2<=161700000 && $12>7.30103 {print $3"\t"$5"\t"$9}' > $W/lpa_all.txt
zcat $M/lipids_LPA.txt.gz | awk 'NR>1 && $1==6 && $2>=160500000 && $2<=161700000 && $12>7.30103 {print $3"\t"$12}' | sort -k2,2gr | awk '{print $1}' > $W/lpa_rank.txt
"/path/to/software/plink" --bfile $SCRATCH/lpa_ldref \
  --clump <(echo -e "SNP\tP"; zcat $M/lipids_LPA.txt.gz | awk 'NR>1 && $1==6 && $2>=160500000 && $2<=161700000 && $12>7.30103 {if($12>m)m=$12; id[NR]=$3; lp[NR]=$12}
      END{lo=7.30103; for(i in id){e=(m<=300)?lp[i]:(lo+(300-lo)*(lp[i]-lo)/(m-lo)); printf "%s\t%.10g\n",id[i],10^-e}}') \
  --clump-p1 5e-8 --clump-p2 1 --clump-r2 0.01 --clump-kb 1000 --out $W/lpa 2>/dev/null || true
awk 'NR>1 && $3!=""{print $3}' $W/lpa.clumped > $W/lpa_keep.txt
grep -wFf $W/lpa_keep.txt $W/lpa_all.txt > $W/lpa_score.txt
echo "  Lp(a) instruments: $(wc -l < $W/lpa_score.txt)"
# LDL-C: genome-wide, distance-pruned at 10 Mb (matches the MVMR polygenic rule)
zcat $M/lipids_LDL_C.txt.gz | awk 'NR>1 && $12>7.30103 {print $1"\t"$2"\t"$3"\t"$5"\t"$9"\t"$12}' | sort -k6,6gr | \
  awk '{k=$1":"int($2/10000000); if(!(k in s)){s[k]=1; print $3"\t"$4"\t"$5}}' > $W/ldl_score.txt
echo "  LDL-C instruments: $(wc -l < $W/ldl_score.txt)"

echo "[2] per-individual scores"
for T in lpa ldl; do
  "${PLINK2}" --bfile "${UKB_GENO_STEP2}" --keep $_REPO/results/agent1_genetics/phenotypes/regenie_keep.txt \
    --extract <(cut -f1 $W/${T}_score.txt) --export A \
    --out $W/dose_$( [ $T = lpa ] && echo LPA || echo LDL ) --threads ${SLURM_NTASKS} --memory 120000
done

echo "[3] statin variable, union of 6153 and 6177"
P=/path/to/ukb_phenotypes
awk -F',' 'FNR==1{next}{gsub(/"/,""); split($0,f,","); e=f[1]; if(f[2]==1)st[e]=1; seen[e]=1}
 END{print "IID\tstatin"; for(x in seen) print x"\t"(x in st?1:0)}' $P/fid6153.csv $P/fid6177.csv > $W/statin.tsv
echo "  statin users: $(awk 'NR>1 && $2==1' $W/statin.tsv | wc -l)"

echo "[4] interaction models"
Rscript $_REPO/scripts/agent1_genetics/R/statin_interaction.R \
  --scores $W/dose_LPA.raw,$W/dose_LDL.raw \
  --weights $W/lpa_score.txt,$W/ldl_score.txt \
  --pheno_q $_REPO/results/agent1_genetics/phenotypes/pheno_quantitative_statinfix.txt \
  --pheno_b $_REPO/results/agent1_genetics/phenotypes/pheno_binary.txt \
  --covar   $_REPO/results/agent1_genetics/phenotypes/covariates.txt \
  --statin  $W/statin.tsv \
  --out     $_REPO/exports/agent1/statin_interaction.csv
echo ALL_DONE
