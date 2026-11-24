#!/usr/bin/env bash
# Draw every figure of revision r4 in the one theme (theme.py): the four main figures, the Supplementary Figure on
# eQTL colocalization, the dose-response Supplementary Figure (written by fig2.py) and the redrawn Supplementary
# Figures in si/. Reads only tables in ../../supplementary, ../../figure_data and ../../../evidence; fits nothing.
set -euo pipefail
cd "$(dirname "$0")"
for s in fig1.py fig2.py fig3.py fig4.py figS_eqtl_coloc.py; do python3 "$s"; done
if [ -d si ]; then for s in si/*.py; do [ -f "$s" ] && python3 "$s"; done; fi
echo "all figures drawn"
