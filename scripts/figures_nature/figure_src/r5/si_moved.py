#!/usr/bin/env python3
"""Two Supplementary Figures (r5) holding r4 main-figure panels that moved to the SI with Supplementary Note 6.

figS_ldl_apob_axis          a  genetic correlation between lipid pairs      (r4 Fig. 1b, unchanged)
                            b  LDL-C and ApoB modelled jointly              (r4 Fig. 1c, unchanged)
figS_ancestry_replication   a  within-ancestry MR against MVP CAD           (r4 Fig. 4b, unchanged)
                            b  cross-ancestry instrument transferability    (r4 Fig. 4c, unchanged)

Reads the deposited tables (source/supplementary) and figure_data/; fits nothing.
"""
import os

import common as C
from common import T
import fig1 as F1   # r4
import fig4 as F4   # r4


def ldl_apob():
    C.relabel("S ldlapob", {("1", "b"): "a", ("1", "c"): "b"})
    P = T.Page(232.0)
    P.panel(2, 10, "a", "LDL-C and ApoB are genetically indistinguishable")
    axa = P.ax(104, 28, 70, 62); F1.draw_b(axa)
    P.panel(262, 10, "b", "Modelled jointly, they cannot be ranked")
    P.key(272, 22, [("s", T.INK, "white", "joint model"), ("s", "white", T.INK, "joint, CI includes 1"), ("o", T.NS, "white", "each fitted alone")], gap=8)
    axb = P.ax(318, 34, 116, 170); F1.draw_c(axb)
    P.save(os.path.join(C.SUP, "figS_ldl_apob_axis"))
    C.write("FigS_ldl_apob_axis_source_data.csv")


def ancestry():
    C.relabel("S ancestryrep", {("4", "b"): "a", ("4", "c"): "b"})
    P = T.Page(214.0)
    P.panel(2, 10, "a", "Within ancestry (MVP CAD)")
    P.key(12, 22, [("s", T.INK, "white", "weighted median"), ("o", T.INK, "white", "IVW"), ("s", "white", T.INK, "open: CI includes 1")], gap=7)
    axb = [P.ax(70, 46, 76, 132), P.ax(158, 46, 76, 132)]; F4.draw_b(axb)
    P.panel(300, 10, "b", "Instrument transferability")
    axc = P.ax(330, 32, 110, 110)
    out = F4.draw_c(axc)
    y = 174
    for anc, mk, r, n in out:
        P.key(334, y, [(mk, "white", T.TRAIT["LDL-C"], "%s, $r$ = %.2f" % (anc, r))], gap=4)
        y += 10
    P.save(os.path.join(C.SUP, "figS_ancestry_replication"))
    C.write("FigS_ancestry_replication_source_data.csv")


if __name__ == "__main__":
    ldl_apob()
    ancestry()
    print("wrote the two moved Supplementary Figures (r5)")
