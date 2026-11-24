#!/bin/bash
#SBATCH -J ukbppp
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 02:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o /path/to/cad-genetics/ukbppp_%j.out
#SBATCH -e /path/to/cad-genetics/ukbppp_%j.err
# ============================================================================
# Cross-platform check on the deCODE cis-pQTL colocalizations.
#
# WHY. The deCODE run (job 3449363) found PCSK9 protein <-> LDL-C <-> CAD all at
# PP.H4 = 1.00 with lead rs11591147 -- a MISSENSE variant. deCODE uses SomaScan
# aptamers, and a coding change can disrupt aptamer binding, mimicking a true
# protein-level signal. Eldjarn et al. Nature 2023 (doi:10.1038/s41586-023-06563-x),
# deCODE's own platform comparison, reports only modest correlation between
# SomaScan and Olink, cis-pQTL support for assay performance in 43% of SomaScan
# versus 72% of Olink assays, and explicit examples where platform choice changes
# disease conclusions. So the deCODE result is not trustworthy until an
# antibody-based platform is checked.
#
# UKB-PPP is Olink Explore 3072 (proximity extension, antibody-based) in ~54,000
# UK Biobank participants. APOC3 is NOT on that panel, so this covers 3 of the 4
# genes: PCSK9, APOB, ANGPTL3.
#
# CAVEAT THIS DOES NOT REMOVE. Olink PEA is also antibody-binding and so is not
# immune to epitope effects; concordance reduces the artifact concern rather than
# eliminating it. And UKB-PPP participants ARE UK Biobank, so unlike deCODE this
# arm shares samples with the outcome GWAS -- acceptable for colocalization, which
# asks about shared causal variants rather than estimating an effect, but it must
# be stated.
#
# SYNAPSE_TOKEN comes from the environment and is never written to disk.
# ============================================================================
set -uo pipefail
_REPO=/path/to/cad-genetics
source $_REPO/config/tacc_paths.sh
set +u; source "${CONDA_BASE}/etc/profile.d/conda.sh"; conda activate "${CONDA_ENV_R}"; set -u
: "${SYNAPSE_TOKEN:?SYNAPSE_TOKEN not set -- resubmit with --export}"
W=$SCRATCH/ukbppp; mkdir -p $W; cd $W
API=https://repo-prod.prod.sagebase.org/repo/v1

declare -A SYN=( [PCSK9]=syn51468818 [APOB]=syn52361466 [ANGPTL3]=syn51468817 )
echo "=== downloading UKB-PPP European pGWAS (Olink) ==="
for G in PCSK9 APOB ANGPTL3; do
  ID=${SYN[$G]}
  [ -s "$G.tar" ] && { echo "  $G.tar present"; continue; }
  curl -sSL --max-time 1800 -H "Authorization: Bearer $SYNAPSE_TOKEN" \
       "$API/entity/$ID/file" -o "$G.tar" || { echo "  $G DOWNLOAD FAILED"; exit 1; }
  echo "  $G.tar $(du -h $G.tar | cut -f1)"
done

echo; echo "=== unpacking and locating cis-window chromosome files ==="
for G in PCSK9 APOB ANGPTL3; do
  mkdir -p $W/$G && tar xf "$G.tar" -C $W/$G 2>/dev/null
  echo "  $G:"; find $W/$G -type f | head -4 | sed 's/^/    /'
  echo "    files: $(find $W/$G -type f | wc -l)"
done
echo ALL_DONE
