#!/usr/bin/env python3
"""Supplementary Figure (figS2_mvmr_pleiotropy): the MVMR-Egger intercept for each outcome, with its 95% CI
(SE derived from the reported P), outcomes ordered by coronary relevance.

Redraw of figure_src/make_figS2.py in the r4 theme. That script listed seven outcomes; Table S3 now holds eight
(peripheral artery disease added), as the figure it replaced did, so all eight are drawn. Reads
source/supplementary/TableS3_mvmr.csv (mvegger_intercept, mvegger_intercept_p); single column (88 mm).
"""
import os
import sys

from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..")))
import theme as T  # noqa: E402

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")

GROUPS = [("Coronary endpoints", [("CAD", "CAD"), ("MI", "Myocardial infarction"), ("REVASC", "Revascularization"), ("MACE", "MACE")]),
          ("Other endpoints", [("PAD", "Peripheral artery disease"), ("STROKE", "Ischaemic stroke"), ("HF", "Heart failure"),
                               ("CV_DEATH", "Cardiovascular death")])]


def intercepts():
    seen = {}
    for r in T.rows(os.path.join(SUP, "TableS3_mvmr.csv")):
        v = (float(r["mvegger_intercept"]), float(r["mvegger_intercept_p"]))
        assert seen.setdefault(r["outcome"], v) == v, r["outcome"]          # one intercept per outcome
    assert sorted(seen) == sorted(o for _, g in GROUPS for o, _ in g), sorted(seen)
    return seen


def build():
    d = intercepts()
    P = T.Page(196.0, w=88 * T.MM)
    P.key(96, 10, [("o", "white", T.INK, "intercept, 95% CI includes 0")], gap=8)
    ax = P.ax(96, 20, 112, 150)
    T.null(ax, 0.0)
    y, ticks, labels, heads = 0.0, [], [], []
    for gi, (head, outs) in enumerate(GROUPS):
        heads.append(y)
        y += 1.0
        for o, lab in outs:
            ic, p = d[o]
            se = abs(ic) / norm.ppf(1 - p / 2) if p < 1 else abs(ic)
            lo, hi = ic - 1.96 * se, ic + 1.96 * se
            T.point(ax, y, ic, lo, hi, T.INK, "o", filled=T.excludes(lo, hi, 0.0))
            ax.text(1.04, y, "P = %.2f" % p, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=T.FS["note"], color=T.MUTED)
            T.src("S2", "", lab, "MVMR-Egger intercept", ic, lo, hi, p="%.3g" % p)
            ticks.append(y); labels.append(lab)
            y += 1.0
        y += 0.3
    for yh, (head, _) in zip(heads, GROUPS):
        ax.text(-0.04, yh + 0.1, head, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=T.FS["note"],
                color=T.MUTED, fontweight="bold")
    ax.set_yticks(ticks); ax.set_yticklabels(labels)
    ax.set_ylim(y - 0.6, -0.5)
    ax.set_xlim(-0.0075, 0.0075); ax.set_xticks([-0.006, -0.003, 0, 0.003, 0.006])
    ax.set_xticklabels(["\u22120.006", "\u22120.003", "0", "0.003", "0.006"])
    ax.set_xlabel("MVMR-Egger intercept (95% CI)")
    T.clean(ax)
    P.save(os.path.join(SRC_DIR, "supplementary", "figS2_mvmr_pleiotropy"))


if __name__ == "__main__":
    build()
    print("wrote supplementary/figS2_mvmr_pleiotropy.pdf")
