#!/usr/bin/env python3
"""Fig. 3 (r5): Lp(a) is independent of the LDL/ApoB axis, and the estimates replicate outside UK Biobank.

a  multivariable MR of LDL-C, triglycerides, HDL-C and Lp(a) against eight endpoints  (r4 Fig. 1d, unchanged)
b  joint MR on CAD in UK Biobank and FinnGen; univariable MR against CARDIoGRAMplusC4D (r4 Fig. 4a, unchanged)

Reads the deposited tables (source/supplementary); fits nothing.
"""
import os

import common as C
from common import T
import fig1 as F1   # r4
import fig4 as F4   # r4


def build():
    C.relabel("3", {("1", "d"): "a", ("4", "a"): "b"})
    P = T.Page(300.0)
    P.panel(2, 10, "a", "Lp(a) is independent of the LDL/ApoB axis")
    P.key(12, 24, [("s", T.TRAIT["LDL-C"], "white", "LDL-C"), ("s", T.TRAIT["TG"], "white", "Triglycerides"),
                   ("s", T.TRAIT["HDL-C"], "white", "HDL-C"), ("s", T.TRAIT["Lp(a)"], "white", "Lp(a)")], gap=8)
    P.key(12, 34, [("s", T.INK, "white", r"$P<1.6\times10^{-3}$"), ("s", T.tint(T.INK, 0.45), "white", r"$P<0.05$"),
                   ("s", "white", T.INK, r"$P\geq0.05$")], gap=8)
    axa = P.ax(104, 46, 150, 222); F1.draw_d(axa)
    P.panel(280, 10, "b", "Replication outside UK Biobank")
    P.key(290, 24, [("s", T.INK, "white", "UK Biobank"), ("^", T.INK, "white", "FinnGen"), ("o", T.INK, "white", "CARDIoGRAMplusC4D")], gap=7)
    P.key(290, 34, [("s", "white", T.INK, "open: 95% CI includes 1")], gap=7)
    axm = P.ax(340, 60, 74, 150); axu = P.ax(426, 60, 74, 150); F4.draw_a(axm, axu)
    P.save(os.path.join(C.SRC_DIR, "figures", "fig3_lpa_axis"))
    C.write("Fig3_source_data.csv")


if __name__ == "__main__":
    build()
    print("wrote Fig. 3 (r5)")
