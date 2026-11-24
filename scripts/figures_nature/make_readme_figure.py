#!/usr/bin/env python3
"""Draw docs/study_design.png, the figure at the top of the README: panel a of the paper's Fig. 1 on its own.

  python scripts/figures_nature/make_readme_figure.py WORK_DIR [OUT_PNG]

WORK_DIR is a folder laid out by make_layout.py: the panel reads the cohort sizes and event counts from
evidence/lpa_emulation/ and the size of the causal-analysis sample from the Supplementary Data. The drawing is done by the
paper's own panel functions, draw_a and draw_a_right in figure_src/r4/fig1.py, unchanged; only the page around them differs.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    work = os.path.abspath(sys.argv[1])
    out = os.path.abspath(sys.argv[2]) if len(sys.argv) > 2 else os.path.join(REPO, "docs", "study_design.png")
    sys.path.insert(0, os.path.join(work, "source", "figure_src", "r4"))
    import theme as T   # noqa: E402  (the copy in WORK_DIR, so that the panel reads WORK_DIR's inputs)
    import fig1 as F1   # noqa: E402
    P = T.Page(126.0)
    F1.draw_a(P, 6, 4)
    F1.draw_a_right(P, 362, 4)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    P.fig.savefig(out, dpi=220)
    print("wrote", out)


if __name__ == "__main__":
    main()
