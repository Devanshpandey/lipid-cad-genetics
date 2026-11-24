#!/usr/bin/env python3
"""Fig. 1 (r5): genetically predicted Lp(a) before and after disease onset, and the events that carry it.

a  study design (r4 Fig. 1a, unchanged)
b  stage contrast for Lp(a), the LDL-C score and the non-lipid score (r4 Fig. 2a, unchanged)
c  Lp(a) on three composites (any, acute-linked or no revascularization) at three stages, beside the trial   (new)
d  Lp(a) by component at three stages, beside the trial's components (r4 Fig. 2b with the third stage and the
   acute-linked revascularization row added)

Reads only evidence/lpa_emulation; fits nothing.
"""
import os

import common as C
from common import T
import fig1 as F1   # r4
import fig2 as F2   # r4

STAGES = (("C_std", "C standardized", "before a first event"), ("A", "A", "established disease"), ("D", "D", "after a first event"))
TRIAL = F2.TRIAL
LP = T.TRAIT["Lp(a)"]


def matched():
    return C.rows(os.path.join(C.EV, "matched_endpoint", "matched_estimates.csv"))


def m_est(M, cohort, endpoint):
    r = C.one(M, cohort=cohort, endpoint=endpoint, horizon="full", model="pooled", exposure="x_lpa")
    return T.num(r["hr_lower"]), T.num(r["lower_lo"]), T.num(r["lower_hi"]), int(r["events"])


# ------------------------------------------------------------------ c: three composites
def groups_c(M):
    out = []
    for head, e in (("Any revascularization", "mace4"), ("Acute-linked only (trial proxy)", "m4u7"), ("No revascularization", "mace3")):
        items = []
        for _, mc, lab in STAGES:
            hr, lo, hi, ev = m_est(M, mc, e)
            items.append((lab, hr, lo, hi, LP, F2.STAGE[{"C standardized": "C_std"}.get(mc, mc)]))
        out.append((head, items))
    out.append(("Lp(a)HORIZON", [("four-point composite", *TRIAL["mace4"], T.TRIAL, "D"),
                                 ("three-point composite", *TRIAL["mace3"], T.TRIAL, "D")]))
    return out


# ------------------------------------------------------------------ d: components
COMPONENTS = [("revasc", "Revascularization, any", None), ("urev7", "Revascularization,\nacute-linked", "revasc"),
              ("mi", "Myocardial infarction", "mi"), ("cvdeath", "Cardiovascular death", "cvdeath"), ("stroke", "Ischaemic stroke", "stroke")]


def draw_d(ax, D, M):
    T.null(ax)
    y, ticks = 0.0, []
    for k, (e, lab, tr) in enumerate(COMPONENTS):
        n = 3 + (tr is not None)
        if k % 2 == 0:
            T.band(ax, y - 0.45, y + (n - 1) * 0.62 + 0.45)
        for j, (cohort, mc, slab) in enumerate(STAGES):
            if e == "urev7":
                hr, lo, hi, ev = m_est(M, mc, e)
            else:
                b, s, _, ev = D["pick"](cohort, e, "x_lpa")
                hr, lo, hi = T.lower(b, s)
            F2.mark(ax, y + j * 0.62, hr, lo, hi, LP, F2.STAGE[cohort])
            T.src("1", "d", lab.replace("\n", " "), slab, hr, lo, hi, events=ev)
        if tr:
            F2.mark(ax, y + 3 * 0.62, *TRIAL[tr], T.TRIAL, "D")
            T.src("1", "d", lab.replace("\n", " "), "Lp(a)HORIZON", *TRIAL[tr])
        ticks.append(y + (n - 1) * 0.62 / 2)
        y += (n - 1) * 0.62 + 1.0
    ax.set_yticks(ticks); ax.set_yticklabels([l for _, l, _ in COMPONENTS], linespacing=1.1)
    ax.set_ylim(y - 1.0 + 0.5, -0.5)
    T.logx(ax, 0.55, 1.45, [0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4])
    ax.set_xlabel("HR per 144 nmol/L lower Lp(a)"); T.clean(ax)


def build():
    C.relabel("1", {("1", "a"): "a", ("2", "a"): "b", ("2", "c"): "c"})
    D, M = F2.load(), matched()
    P = T.Page(476.0)
    P.panel(2, 10, "a", "Study design")
    F1.draw_a(P, 6, 16)
    F1.draw_a_right(P, 362, 16)
    F2.stage_key(P, 10, 142)
    P.key(10, 153, [("line", None, T.TRAIT["Lp(a)"], "Lp(a)"), ("line", None, T.TRAIT["LDL-C"], "LDL-C score (positive control)"),
                    ("line", None, T.TRAIT["non-lipid"], "non-lipid coronary score")], gap=9)
    P.panel(2, 170, "b", "The associations persist after disease onset")
    axb = P.ax(80, 184, 72, 104); F2.draw_a(axb, D)
    P.panel(2, 322, "c", "How much persists depends on the events counted")
    axc = P.ax(118, 336, 72, 118)
    F2.grouped(axc, groups_c(M), "c", (0.66, 1.22), [0.7, 0.8, 0.9, 1, 1.2], "HR per 144 nmol/L lower Lp(a)")
    P.panel(262, 170, "d", "The persisting association is carried by coronary events")
    axd = P.ax(352, 184, 152, 270); draw_d(axd, D, M)
    P.save(os.path.join(C.SRC_DIR, "figures", "fig1_lpa_stages"))
    C.write("Fig1_source_data.csv")


if __name__ == "__main__":
    build()
    print("wrote Fig. 1 (r5)")
