#!/usr/bin/env python3
"""Supplementary Figure (moved from r3 Fig. 1b): LDL-C and ApoB make the same colocalization call with gene expression.

One point per gene x tissue cis-eQTL test: posterior of a shared causal variant (PP.H4) with LDL-C against that with ApoB.
Reads source/supplementary/TableS27b_eqtl_coloc_ldl_apob.csv; single column (88 mm).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import theme as T  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.abspath(os.path.join(HERE, "..", ".."))


def build():
    d = {}
    for r in T.rows(os.path.join(SRC_DIR, "supplementary", "TableS27b_eqtl_coloc_ldl_apob.csv")):
        d.setdefault((r["gene"], r["tissue"]), {})[r["trait"]] = float(r["PP.H4"])
    pairs = [(k, v["lipids_LDL_C"], v["lipids_APOB"]) for k, v in d.items() if "lipids_LDL_C" in v and "lipids_APOB" in v]
    disc = [p for p in pairs if (p[1] > 0.8) != (p[2] > 0.8)]
    conc = [p for p in pairs if (p[1] > 0.8) == (p[2] > 0.8)]
    assert len(pairs) == 78 and len(disc) == 1, (len(pairs), len(disc))
    P = T.Page(200.0, w=88 * T.MM)
    ax = P.ax(34, 14, 200, 150)
    ax.plot([0, 1], [0, 1], color=T.NULL, lw=T.LW["ref"], zorder=1)
    ax.axhline(0.8, color=T.BAND, lw=1.2, zorder=0); ax.axvline(0.8, color=T.BAND, lw=1.2, zorder=0)
    ax.plot([p[1] for p in conc], [p[2] for p in conc], ls="none", marker="o", ms=3.0, mfc=T.TRAIT["LDL-C"], mec="white", mew=0.4, alpha=0.85, zorder=3)
    for (g, tis), x, y in disc:
        ax.plot([x], [y], ls="none", marker="o", ms=3.6, mfc="white", mec=T.INK, mew=0.9, zorder=4)
        ax.annotate("$\\it{%s}$, %s" % (g, tis.replace("Artery_", "").replace("_", " ").lower()), (x, y), xytext=(0.52, 0.12),
                    fontsize=T.FS["note"], arrowprops=dict(arrowstyle="-", lw=0.4, color=T.MUTED, shrinkA=0.5, shrinkB=2))
        T.src("S", "eqtl", "%s %s" % (g, tis), "LDL-C vs ApoB PP.H4", x, y, y)
    ax.text(0.02, 0.98, "%d tests; same call at 0.8 in %d" % (len(pairs), len(conc)), transform=ax.transAxes, ha="left", va="top", fontsize=T.FS["note"], color=T.MUTED)
    ax.set_xlim(-0.04, 1.04); ax.set_ylim(-0.04, 1.04)
    ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1]); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_xlabel("PP.H4 with LDL-C"); ax.set_ylabel("PP.H4 with ApoB")
    T.clean(ax, left=True)
    P.save(os.path.join(SRC_DIR, "supplementary", "figS_eqtl_coloc"))


if __name__ == "__main__":
    build()
    print("wrote the eQTL colocalization Supplementary Figure")
