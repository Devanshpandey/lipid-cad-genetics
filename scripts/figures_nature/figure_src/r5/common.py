"""Shared set-up for the r5 figures: the r4 theme and panel functions, unchanged, with source data relabelled.

r5 reorders the main figures (the comparison with Lp(a)HORIZON first) and moves two r4 panels to the Supplementary
Information. The theme (r4/theme.py) and the r4 panel functions are imported, not copied, so every figure keeps one
style. The r4 functions tag their source-data rows with r4 figure numbers; relabel() maps them to the r5 ones.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
R4 = os.path.join(os.path.dirname(HERE), "r4")
sys.path.insert(0, R4)
import theme as T  # noqa: E402

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")
EV = os.path.abspath(os.path.join(SRC_DIR, "..", "evidence", "lpa_emulation"))
OUT_DATA = os.path.join(HERE, "source_data")

_orig_src = T.src


def relabel(figure, panels=None):
    """Send every T.src call to figure `figure`, mapping r4 (figure, panel) pairs to r5 panel letters."""
    panels = panels or {}
    T.SRC.clear()

    def src(fig, panel, *a, **k):
        _orig_src(figure, panels.get((fig, panel), panel), *a, **k)
    T.src = src


def write(name):
    os.makedirs(OUT_DATA, exist_ok=True)
    T.write_src(os.path.join(OUT_DATA, name))


def one(table, **kw):
    r = [x for x in table if all(x[k] == v for k, v in kw.items())]
    assert len(r) == 1, (kw, len(r))
    return r[0]


def rows(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))
