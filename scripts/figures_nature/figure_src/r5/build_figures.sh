#!/usr/bin/env bash
# Draw every figure of revision r5 in the one theme (../r4/theme.py).
#   Main figures: r5/fig1.py, r5/fig2.py, r5/fig3.py and r4/fig3.py (Fig. 4, unchanged from r4 Fig. 3).
#   Supplementary Figures: r5/si_moved.py (two panels pairs moved from r4 Figs. 1 and 4), r4/fig2.py (the dose-response
#   Supplementary Figure; it also redraws the unused r4 Fig. 2), r4/figS_eqtl_coloc.py and r4/si/*.py.
# Reads only tables in ../../supplementary, ../../figure_data and ../../../evidence; fits nothing.
set -euo pipefail
cd "$(dirname "$0")"
( cd ../r4 && for s in fig2.py fig3.py figS_eqtl_coloc.py; do python3 "$s"; done
  if [ -d si ]; then for s in si/*.py; do [ -f "$s" ] && python3 "$s"; done; fi )
for s in fig1.py fig2.py fig3.py si_moved.py; do python3 "$s"; done
echo "all r5 figures drawn"
