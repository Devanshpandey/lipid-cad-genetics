#!/usr/bin/env python3
"""Lay out the inputs of the figure code as the manuscript folder had them, and optionally draw the figures.

  python scripts/figures_nature/make_layout.py WORK_DIR --supplementary-data DIR [--figure-data DIR] [--draw]

The figure code in figure_src/ was run inside the manuscript's folder and reads its inputs by relative path:
  source/supplementary/      the paper's Supplementary Data files, under the working names the code uses
  source/figure_data/        a few further summary tables (one is in this repository; see INPUTS.tsv)
  evidence/lpa_emulation/    the comparison with Lp(a)HORIZON (in this repository: results_summary/lpa_horizon/)
This script copies the code and those inputs into WORK_DIR in that layout, renaming each Supplementary Data file to
its working name and checking every input's md5 against INPUTS.tsv. With --draw it runs the figure scripts in the
order of figure_src/r5/build_figures.sh and reports any it could not run. Figures land in WORK_DIR/source/figures
(main) and WORK_DIR/source/supplementary (Supplementary). Nothing is written inside the repository.

DIR for --supplementary-data is a folder holding the paper's Supplementary Data files as published
(Supplementary_Table_*.csv and Supplementary_Data_*.csv). --figure-data is only needed for the Supplementary Figures
that read the authors' figure-input copies listed in INPUTS.tsv as "not in this repository".
"""
import argparse
import csv
import hashlib
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
ORDER = [("r4", "fig2.py"), ("r4", "fig3.py"), ("r4", "figS_eqtl_coloc.py"), ("r4", "si"),
         ("r5", "fig1.py"), ("r5", "fig2.py"), ("r5", "fig3.py"), ("r5", "si_moved.py")]


def md5(p):
    with open(p, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("work")
    ap.add_argument("--supplementary-data", required=True)
    ap.add_argument("--figure-data")
    ap.add_argument("--draw", action="store_true")
    a = ap.parse_args()
    W = os.path.abspath(a.work)
    src = os.path.join(W, "source")
    shutil.copytree(os.path.join(HERE, "figure_src"), os.path.join(src, "figure_src"), dirs_exist_ok=True)
    for d in ("supplementary", "figure_data", "figures"):
        os.makedirs(os.path.join(src, d), exist_ok=True)

    missing, wrong, n = [], [], 0
    with open(os.path.join(HERE, "INPUTS.tsv"), newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["origin"] == "repository":
                s = os.path.join(REPO, r["path_in_origin"])
            elif r["origin"] == "Supplementary Data of the paper":
                s = os.path.join(a.supplementary_data, r["path_in_origin"])
            else:
                s = os.path.join(a.figure_data, r["path_in_origin"]) if a.figure_data else None
            if not s or not os.path.exists(s):
                missing.append(r["working_path"])
                continue
            dst = os.path.join(W, r["working_path"])
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copyfile(s, dst)
            n += 1
            if md5(dst) != r["md5"]:
                wrong.append(r["working_path"])
    print("inputs laid out: %d; missing: %d; md5 differs from the paper's: %d" % (n, len(missing), len(wrong)))
    for x in missing:
        print("  missing:", x)
    for x in wrong:
        print("  md5 differs:", x)
    if not a.draw:
        return 0

    failed = []
    for sub, s in ORDER:
        d = os.path.join(src, "figure_src", sub)
        scripts = [os.path.join("si", f) for f in sorted(os.listdir(os.path.join(d, "si"))) if f.endswith(".py")] if s == "si" else [s]
        for sc in scripts:
            r = subprocess.run([sys.executable, sc], cwd=d, capture_output=True, text=True)
            tag = "%s/%s" % (sub, sc)
            if r.returncode:
                failed.append(tag)
                print("  not drawn: %s (%s)" % (tag, (r.stderr.strip().splitlines() or ["?"])[-1][:120]))
            else:
                print("  drawn:", tag)
    print("scripts run: %d failed" % len(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
