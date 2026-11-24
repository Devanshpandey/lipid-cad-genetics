#!/usr/bin/env python3
"""Supplementary Figure (PDE3B rare coding variation), redrawn in the r4 theme.

a  phenome-wide scan of 1,251 pLoF carriers: quantitative phenotypes (SD) and binary outcomes (OR), 95% CI;
   filled = Bonferroni-significant
b  triglyceride effect by variant class (pLoF > missense > synonymous negative control)
c  CAD association by design: gene-burden proxies (OR) and the prospective carrier Cox model (HR)

Same content as figure_src/make_fig4_pde3b.py. Reads source/supplementary/TableS13_pde3b_phewas.csv,
pde3b_consequence_class.csv and TableS20_pde3b_cox_cad.csv. The two gene-burden proxy estimates in c were typed into the
old script and are not in any deposited table; they are kept as typed. The Cox estimate is read from Table S20 and checked
against the value the old script typed.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import theme as T  # noqa: E402

SRC_DIR = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SUP = os.path.join(SRC_DIR, "supplementary")
TG = T.TRAIT["TG"]
AXIS = {"TRIGLY": T.TRAIT["TG"], "HDL_C": T.TRAIT["HDL-C"], "APOA1": T.TRAIT["HDL-C"], "APOB": T.TRAIT["ApoB"],
        "nonHDL_C": T.TRAIT["LDL-C"], "LDL_C": T.TRAIT["LDL-C"]}
NICE = {"TRIGLY": "Triglycerides", "HDL_C": "HDL-C", "APOA1": "ApoA1", "APOB": "ApoB", "nonHDL_C": "non-HDL-C", "LDL_C": "LDL-C",
        "CRP": "CRP", "eGFR": "eGFR", "HBA1C": "HbA1c", "STATIN_USE": "Statin use", "CAD": "CAD", "MI": "MI", "MACE": "MACE",
        "REVASC": "Revascularization", "PAD": "PAD", "STROKE": "Stroke", "HF": "Heart failure", "AF": "Atrial fibrillation",
        "CV_DEATH": "CV death"}
QUANT = ["TRIGLY", "HDL_C", "APOA1", "APOB", "nonHDL_C", "LDL_C", "CRP", "eGFR", "HBA1C"]
BINARY = ["STATIN_USE", "CAD", "MI", "MACE", "REVASC", "PAD", "STROKE", "HF", "AF", "CV_DEATH"]


def ci(r):
    eff, se = float(r["effect"]), float(r["se"])
    if r["effect_lab"] == "OR":
        return eff, math.exp(math.log(eff) - 1.96 * se), math.exp(math.log(eff) + 1.96 * se)
    return eff, eff - 1.96 * se, eff + 1.96 * se


def rowsetup(ax, labels):
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels)
    ax.set_ylim(len(labels) - 0.5, -0.5)
    T.clean(ax)


def draw_a(axq, axb):
    ph = {r["pheno"]: r for r in T.rows(os.path.join(SUP, "TableS13_pde3b_phewas.csv"))}
    out = []
    T.null(axq, 0.0)
    for i, p in enumerate(QUANT):
        e, lo, hi = ci(ph[p]); sig = float(ph[p]["P_bonf"]) < 0.05
        T.point(axq, i, e, lo, hi, AXIS.get(p, T.TRAIT["non-lipid"]), "o", filled=sig)
        out.append((p, round(e, 3), round(lo, 3), round(hi, 3), sig))
    rowsetup(axq, [NICE[p] for p in QUANT])
    axq.set_xlim(-0.48, 0.42); axq.set_xticks([-0.4, -0.2, 0, 0.2, 0.4])
    axq.set_xticklabels(["−0.4", "−0.2", "0", "0.2", "0.4"])
    axq.set_xlabel("Effect (SD, 95% CI)")
    T.null(axb, 1.0)
    for i, p in enumerate(BINARY):
        e, lo, hi = ci(ph[p]); sig = float(ph[p]["P_bonf"]) < 0.05
        T.point(axb, i, e, lo, hi, T.TRAIT["non-lipid"] if p == "STATIN_USE" else T.INK, "o", filled=sig)
        out.append((p, round(e, 3), round(lo, 3), round(hi, 3), sig))
    rowsetup(axb, [NICE[p] for p in BINARY])
    T.logx(axb, 0.5, 1.8, [0.5, 0.7, 1, 1.4])
    axb.set_xlabel("Odds ratio (95% CI)")
    return out


def draw_b(ax):
    cc = T.rows(os.path.join(SUP, "pde3b_consequence_class.csv"))
    labs = {"pLoF": "pLoF", "missense": "missense", "synonymous": "synonymous\n(neg. control)"}
    T.null(ax, 0.0)
    out = []
    for i, r in enumerate(cc):
        b, lo, hi = float(r["beta"]), float(r["lo"]), float(r["hi"])
        T.point(ax, i, b, lo, hi, TG, "o", filled=T.excludes(lo, hi, 0.0))
        ax.text(1.03, i, "n = %s" % format(int(r["n_carriers"]), ","), transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=T.FS["note"], color=T.MUTED)
        out.append((r["class"], round(b, 3), round(lo, 3), round(hi, 3), r["n_carriers"]))
    ax.set_yticks(range(len(cc))); ax.set_yticklabels([labs[r["class"]] for r in cc], linespacing=1.0)
    ax.set_ylim(len(cc) - 0.5, -0.5)
    ax.set_xlim(-0.47, 0.16); ax.set_xticks([-0.4, -0.3, -0.2, -0.1, 0, 0.1])
    ax.set_xticklabels(["−0.4", "−0.3", "−0.2", "−0.1", "0", "0.1"])
    ax.set_xlabel("Triglyceride effect (SD, 95% CI)")
    T.clean(ax)
    return out


def draw_c(ax):
    cox = T.rows(os.path.join(SUP, "TableS20_pde3b_cox_cad.csv"))
    assert len(cox) == 1
    hr = tuple(float(cox[0][k]) for k in ("HR", "CI_lo95", "CI_hi95"))
    assert hr == (1.028, 0.838, 1.261) and cox[0]["carriers"] == "1103" and cox[0]["carrier_CAD_events"] == "92", (hr, cox[0])
    designs = [("pLoF-burden proxy", (0.861, 0.713, 1.038), "OR", "o"),
               ("broad-mask proxy", (0.844, 0.704, 1.012), "OR", "o"),
               ("carrier incident CAD (Cox)", hr, "HR", "s")]
    ax.axvspan(0.75, 1.25, color=T.BAND, lw=0, zorder=0)
    T.null(ax, 1.0)
    for i, (nm, (e, lo, hi), lab, mk) in enumerate(designs):
        T.point(ax, i, e, lo, hi, T.INK, mk, filled=T.excludes(lo, hi))
        ax.text(1.03, i, "%s %s" % (lab, T.ci_text(e, lo, hi)), transform=ax.get_yaxis_transform(), ha="left", va="center",
                fontsize=T.FS["note"])
    ax.set_yticks(range(len(designs)))
    ax.set_yticklabels([d[0].replace("incident CAD (Cox)", "incident\nCAD (Cox)") for d in designs], linespacing=1.0)
    ax.set_ylim(len(designs) - 0.5, -0.5)
    T.logx(ax, 0.6, 1.55, [0.6, 0.8, 1, 1.25, 1.5])
    ax.set_xlabel("Effect on CAD (95% CI)")
    T.clean(ax)
    return [(d[0], d[1]) for d in designs]


def build():
    P = T.Page(244.0)
    P.panel(2, 10, "a", "Quantitative phenotypes")
    P.text(11, 19.5, "1,251 pLoF carriers · filled = Bonferroni-significant", color=T.MUTED)
    P.text(154, 10, "Binary outcomes", fontsize=T.FS["title"], fontweight="bold")
    P.text(154, 19.5, "no adverse outcome passes the multiple-", color=T.MUTED)
    P.text(154, 27.0, "testing threshold; CIs remain wide", color=T.MUTED)
    axq = P.ax(46, 42, 86, 172); axb = P.ax(204, 42, 66, 172)
    a_out = draw_a(axq, axb)
    P.panel(290, 10, "b", "Effect scales with variant severity")
    P.text(299, 19.5, "carrier sets analysed separately; synonymous = negative control", color=T.MUTED)
    axs = P.ax(352, 36, 118, 58); b_out = draw_b(axs)
    P.panel(290, 126, "c", "CAD protection not established")
    P.text(299, 135.5, "burden proxy protective but underpowered; carrier Cox null", color=T.MUTED)
    P.key(299, 148, [("o", "white", T.INK, "OR, gene-burden proxy"), ("s", "white", T.INK, "HR, carrier Cox"),
                     ("s", "white", T.INK, "open: CI includes 1")], gap=7)
    P.text(299, 157, "92 events / 1,103 carriers · shaded = effect-ratio band 0.75–1.25 (OR/HR)", color=T.MUTED)
    axc = P.ax(352, 166, 94, 50); c_out = draw_c(axc)
    P.save(os.path.join(SRC_DIR, "figures", "fig4_pde3b"))
    print("a:", a_out); print("b:", b_out); print("c:", c_out)


if __name__ == "__main__":
    build()
