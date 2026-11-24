#!/bin/bash
#------------------------------------------------------------
# Fetch UKB-PPP (Olink Explore 3072) European-discovery pGWAS for the lipid
# effector genes, to pair with the deCODE SomaScan arm.
#
# CREDENTIAL: submit with
#   sbatch --export=ALL,SYNAPSE_TOKEN=<token> 08ab_ukbppp_fetch.slurm
# Never written to disk. Entity IDs were resolved from the European (discovery)
# folder syn51365303 (2,940 files).
#
# Coverage note: APOC3, CETP, LIPC, APOA5, HMGCR, NPC1L1 and TRIB1 are NOT on
# the Olink panel; CETP/LIPC/NPC1L1/TRIB1 are on neither platform, so those
# genes get no protein layer at all.
#
# UKB-PPP participants ARE UK Biobank, so this arm shares samples with the
# outcome GWAS -- acceptable for colocalization (which asks about a shared
# causal variant, not an effect size) but it must be stated.
#------------------------------------------------------------
#SBATCH -J ukbppp_fetch
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t 06:00:00
#SBATCH -A YOUR_ALLOCATION
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err

set -uo pipefail
: "${SYNAPSE_TOKEN:?SYNAPSE_TOKEN not set -- submit with --export=ALL,SYNAPSE_TOKEN=...}"
API=https://repo-prod.prod.sagebase.org/repo/v1
W="${SCRATCH}/pqtl2/ukbppp"; mkdir -p "$W"; cd "$W"

read -r -d '' ENTS <<'E' || true
LDLR syn51469717
LPA syn52361776
SORT1 syn51469328
APOE syn52362532
APOC1 syn52361316
LPL syn51469090
LDLRAP1 syn52361838
APOA1 syn52362862
LCAT syn52363655
ANGPTL4 syn51468799
PLTP syn51469473
E

echo "$ENTS" | while read -r G ID; do
  [ -z "$G" ] && continue
  if [ -s "$G.tar" ]; then echo "  $G.tar present ($(du -h "$G.tar" | cut -f1))"; continue; fi
  curl -sSL --max-time 3600 -H "Authorization: Bearer $SYNAPSE_TOKEN" \
       "$API/entity/$ID/file" -o "$G.tar" \
    && echo "  $G.tar $(du -h "$G.tar" | cut -f1)" \
    || echo "  $G DOWNLOAD FAILED"
done

echo; echo "=== unpacking ==="
for T in *.tar; do
  G="${T%.tar}"; mkdir -p "$G" && tar xf "$T" -C "$G" 2>/dev/null
  echo "  $G: $(find "$G" -type f | wc -l) files"
done
echo ALL_DONE
