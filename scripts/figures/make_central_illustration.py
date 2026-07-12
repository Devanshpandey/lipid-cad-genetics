#!/usr/bin/env python3
"""
make_central_illustration.py — evidence-funnel Central Illustration.

Left-to-right funnel in four stages that mirror the paper's logic:
  1. Discovery data (UK Biobank)                -> four lipid exposure axes
  2. Genetic evidence, three tiers              -> trait / locus / gene
  3. Replication and generalization (full col)  -> EU cohorts, GLGC, MVP ancestries
  4. Gene and candidate prioritization          -> established set + candidates

Design grammar (figstyle): fixed biological colours, readable type. No plaque
drawing, no hexagon, no drug/therapy row. Every number is real.

Output: fig_central_illustration.{png,pdf}
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from figstyle import INK, LIPID, SUBTLE

W, H = 182, 104
fig, ax = plt.subplots(figsize=(13.6, 7.77))
ax.set_xlim(0, W); ax.set_ylim(0, H); ax.set_aspect("equal"); ax.axis("off")

RED, ORNG, BLUE, PURP = LIPID["LDL-C"], LIPID["TG"], LIPID["HDL-C"], LIPID["Lp(a)"]
TEAL, BURG, NS = LIPID["rare"], LIPID["outcome"], LIPID["ns"]
NAVY = "#173A5E"; MID = SUBTLE; LINE = "#D3DAE2"
GOOD = "#2A9D8F"          # supported claim (teal-green check)


def draw_check(x, y, s=1.0, color=GOOD):
    ax.plot([x - 1.1 * s, x - 0.25 * s, x + 1.25 * s],
            [y + 0.0 * s, y - 0.95 * s, y + 1.05 * s],
            color=color, lw=1.7, solid_capstyle="round",
            solid_joinstyle="round", zorder=6)


def rrect(x, y, w, h, fc, ec="none", lw=0, r=2.0, z=2):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle=f"round,pad=0,rounding_size={r}",
                 fc=fc, ec=ec, lw=lw, zorder=z, mutation_aspect=1))


def stage(x, w, header, hcol, tint):
    rrect(x, 9, w, 79, tint, ec=hcol, lw=1.1, r=2.8, z=1)
    rrect(x + 0.7, 82.4, w - 1.4, 5.0, hcol, r=2.0, z=2)
    ax.text(x + w / 2, 84.9, header, ha="center", va="center",
            fontsize=10.4, fontweight="bold", color="white", zorder=3)


def arrow(x1, x2, y=48.5):
    ax.add_patch(FancyArrowPatch((x1, y), (x2, y), arrowstyle="-|>",
                 mutation_scale=20, lw=2.6, color=INK, zorder=7,
                 shrinkA=0, shrinkB=0, capstyle="round"))


# ===================================================================== header
ax.text(W / 2, 100.3,
        "From causal lipid fractions to replicated loci, effector genes and candidate targets",
        ha="center", va="center", fontsize=15.8, fontweight="bold", color=INK)
ax.text(W / 2, 95.0,
        "Common- and rare-variant human genetics of the lipid–coronary artery disease axis, "
        "tested out of sample across cohorts and ancestries",
        ha="center", va="center", fontsize=9.4, color=MID)

X1, W1 = 2.0, 39.0
X2, W2 = 44.5, 42.0
X3, W3 = 90.0, 43.0
X4, W4 = 137.0, 43.0
arrow(X1 + W1, X2)
arrow(X2 + W2, X3)
arrow(X3 + W3, X4)

# ============================================================ 1 · discovery data
stage(X1, W1, "1 · DISCOVERY", NAVY, "#EEF2F8")
cx1 = X1 + W1 / 2
rrect(X1 + 2.4, 71.5, W1 - 4.8, 9.6, "white", ec="#C3CDDA", lw=1.0, r=1.6, z=3)
ax.text(cx1, 79.2, "UK Biobank", ha="center", fontsize=8.9, fontweight="bold", color=INK)
ax.text(cx1, 76.3, "402,200 genotyped participants", ha="center", fontsize=7.3, color=INK)
ax.text(cx1, 73.9, "469,835 sequenced exomes", ha="center", fontsize=7.3, color=INK)
ax.text(cx1, 68.4, "11 biomarkers  ·  10 cardiovascular outcomes",
        ha="center", fontsize=6.9, color=MID)

ax.text(cx1, 63.2, "Lipid exposure axes", ha="center", fontsize=7.6,
        fontweight="bold", color=INK)
axes4 = [
    ("LDL-C / ApoB", RED, 57.0),
    ("Triglyceride-rich / remnant", ORNG, 49.0),
    ("Lipoprotein(a)", PURP, 41.0),
    ("HDL-C / ApoA1", BLUE, 33.0),
]
for name, col, yy in axes4:
    rrect(X1 + 3.0, yy - 2.7, W1 - 6.0, 5.4, "white", ec=col, lw=1.5, r=1.4, z=4)
    ax.add_patch(FancyBboxPatch((X1 + 3.6, yy - 2.1), 1.5, 4.2,
                 boxstyle="round,pad=0,rounding_size=0.5", fc=col, ec="none", zorder=5))
    ax.text(X1 + 6.4, yy, name, ha="left", va="center", fontsize=7.6,
            fontweight="bold", color=INK, zorder=6)
ax.text(cx1, 25.5, "biomarker/renal traits retained\nas controls (not shown)",
        ha="center", va="center", fontsize=6.2, color=MID, fontstyle="italic",
        linespacing=1.15)

# ==================================================== 2 · genetic evidence tiers
stage(X2, W2, "2 · GENETIC EVIDENCE", "#3F6098", "#EDF1F7")
tiers = [
    ("TRAIT-LEVEL", "genetic correlation · multivariable MR · MR-BMA",
     "LDL-C independent causal effect on CAD (OR 1.30/SD);\nLDL-C and ApoB genetically inseparable", 74.5),
    ("LOCUS-LEVEL", "colocalization · SuSiE fine-mapping",
     "five shared lipid–CAD loci (PP.H4 > 0.8):\nPCSK9, SORT1, APOE, LDLR, TRIB1", 53.0),
    ("GENE-LEVEL", "rare-variant burden · allelic-direction",
     "canonical effector genes recovered; LoF directions\nmatch approved lipid therapeutics", 31.5),
]
for label, methods, out, ty in tiers:
    rrect(X2 + 2.2, ty - 12.0, W2 - 4.4, 16.6, "white", ec="#C7D2E2", lw=1.0, r=1.8, z=3)
    rrect(X2 + 2.2, ty + 1.3, W2 - 4.4, 3.3, "#DCE6F3", r=1.0, z=4)
    ax.text(X2 + 4.0, ty + 2.95, label, ha="left", va="center", fontsize=8.0,
            fontweight="bold", color=NAVY, zorder=5)
    ax.text(X2 + 4.0, ty - 1.4, methods, ha="left", va="center", fontsize=7.0,
            color="#3F6098", fontstyle="italic", zorder=5)
    ax.text(X2 + 4.0, ty - 7.4, out, ha="left", va="center", fontsize=7.0,
            color=INK, linespacing=1.25, zorder=5)

# ================================================ 3 · replication & generalization
stage(X3, W3, "3 · REPLICATION", TEAL, "#E9F5F1")
reps = [
    ("European replication", "FinnGen  +  CARDIoGRAMplusC4D", 77.0),
    ("Cross-ancestry transferability", "Global Lipids Genetics Consortium", 70.0),
    ("Within-ancestry causal MR", "African & Hispanic/admixed  ·  MVP CAD", 63.0),
]
for title, detail, yy in reps:
    ax.add_patch(FancyBboxPatch((X3 + 2.6, yy - 1.0), 1.4, 4.6,
                 boxstyle="round,pad=0,rounding_size=0.5", fc=TEAL, ec="none", zorder=4))
    ax.text(X3 + 5.2, yy + 2.4, title, ha="left", fontsize=7.7, fontweight="bold", color=INK)
    ax.text(X3 + 5.2, yy - 0.1, detail, ha="left", fontsize=6.9, color=MID)

# checkmark table: conclusion x (European rep, cross-ancestry)
ty0 = 52.5
ax.text(X3 + W3 / 2, 56.6, "Does the conclusion hold?", ha="center",
        fontsize=7.4, fontweight="bold", color=INK)
colx = [X3 + 21.5, X3 + 31.5, X3 + 39.5]  # label / EU / cross-ancestry
ax.text(colx[1], ty0 + 1.6, "EU", ha="center", fontsize=6.6, fontweight="bold", color=MID)
ax.text(colx[2], ty0 + 1.6, "cross-anc.", ha="center", fontsize=6.6, fontweight="bold", color=MID)
rowdat = [
    ("LDL/ApoB axis", RED, "check", "check"),
    ("Triglycerides", ORNG, "check", "check"),
    ("Lp(a)", PURP, "check", "n.e."),
    ("HDL-C target-independent", BLUE, "incons", "incons"),
]
ry = ty0 - 1.8
for name, col, eu, ca in rowdat:
    ax.add_patch(FancyBboxPatch((X3 + 2.6, ry - 0.9), 1.1, 2.4,
                 boxstyle="round,pad=0,rounding_size=0.4", fc=col, ec="none", zorder=4))
    ax.text(X3 + 4.6, ry, name, ha="left", va="center", fontsize=6.9, color=INK)
    for cxx, val in ((colx[1], eu), (colx[2], ca)):
        if val == "check":
            draw_check(cxx, ry, s=1.0)
        elif val == "n.e.":
            ax.text(cxx, ry, "not\nevaluated", ha="center", va="center", fontsize=5.0,
                    color=MID, fontstyle="italic", linespacing=0.95)
        else:
            ax.text(cxx, ry, "incons.", ha="center", va="center", fontsize=6.0,
                    color=MID, fontstyle="italic")
    ry -= 3.9
draw_check(X3 + 8.2, 33.6, s=0.8)
ax.text(X3 + 9.4, 33.6, "check = supported   ·   incons. = inconsistent across models/ancestry",
        ha="left", va="center", fontsize=5.9, color=MID)
rrect(X3 + 2.6, 12.5, W3 - 5.2, 17.5, "white", ec=TEAL, lw=1.1, r=1.8, z=3)
ax.text(X3 + W3 / 2, 26.8, "Causal architecture reproduces", ha="center",
        fontsize=7.5, fontweight="bold", color=INK)
ax.text(X3 + W3 / 2, 22.6, "LDL-C median OR 1.39 (African)\nOR 1.63 (Hispanic/admixed) MVP CAD",
        ha="center", va="center", fontsize=6.8, color=INK, linespacing=1.25)
ax.text(X3 + W3 / 2, 16.0, "triglyceride effect concordant; HDL-C null",
        ha="center", va="center", fontsize=6.6, color=MID, fontstyle="italic")

# =============================================== 4 · gene / candidate prioritization
stage(X4, W4, "4 · PRIORITIZATION", BURG, "#F6EBEE")
ax.text(X4 + W4 / 2, 78.6, "Established architecture", ha="center",
        fontsize=7.8, fontweight="bold", color=INK)
established = ["PCSK9", "LDLR", "SORT1", "APOE", "APOB", "APOC3", "ANGPTL3", "LPA"]
gx0, gw, gh = X4 + 3.4, 8.7, 4.4
for i, g in enumerate(established):
    col, row = i % 4, i // 4
    gx = gx0 + col * (gw + 0.9)
    gy = 72.0 - row * (gh + 1.2)
    rrect(gx, gy, gw, gh, "white", ec=BURG, lw=1.2, r=1.2, z=4)
    ax.text(gx + gw / 2, gy + gh / 2, g, ha="center", va="center",
            fontsize=6.9, fontweight="bold", color=INK, zorder=5, fontstyle="italic")

rrect(X4 + 2.6, 46.5, W4 - 5.2, 9.6, "#FBF3E4", ec="#C79432", lw=1.1, r=1.8, z=3)
ax.text(X4 + W4 / 2, 53.7, "Integrated common + rare score", ha="center",
        fontsize=7.5, fontweight="bold", color="#9A6E12")
ax.text(X4 + W4 / 2, 49.4,
        "improves genome-wide recovery of established\ngenes: AUROC 0.97 vs 0.90 (burden alone)",
        ha="center", va="center", fontsize=6.8, color=INK, linespacing=1.2)

ax.text(X4 + W4 / 2, 41.4, "Candidate prioritization", ha="center",
        fontsize=7.8, fontweight="bold", color=INK)
cands = ["PDE3B", "ANGPTL8", "PLIN1"]
cw = 11.2
cx0 = X4 + (W4 - (len(cands) * cw + (len(cands) - 1) * 1.4)) / 2
for i, g in enumerate(cands):
    gx = cx0 + i * (cw + 1.4)
    rrect(gx, 35.2, cw, 4.6, "white", ec="#9A6E12", lw=1.3, r=1.2, z=4)
    ax.text(gx + cw / 2, 37.5, g, ha="center", va="center", fontsize=7.1,
            fontweight="bold", color="#9A6E12", zorder=5, fontstyle="italic")

rrect(X4 + 2.6, 13.0, W4 - 5.2, 17.8, "white", ec=BURG, lw=1.1, r=1.8, z=3)
ax.text(X4 + W4 / 2, 27.0, "PDE3B", ha="center", fontsize=8.0, fontweight="bold",
        color=BURG, fontstyle="italic")
ax.text(X4 + W4 / 2, 22.6, "Strong lipid-modifying evidence", ha="center",
        fontsize=7.0, color=INK)
ax.text(X4 + W4 / 2, 18.4, "CAD protection not established", ha="center",
        fontsize=7.2, fontweight="bold", color=BURG)
ax.text(X4 + W4 / 2, 14.8, "(no cis-QTL; carrier incident-CAD null, underpowered)",
        ha="center", fontsize=6.0, color=MID, fontstyle="italic")

fig.subplots_adjust(left=0.004, right=0.996, top=0.996, bottom=0.004)
for ext in ("png", "pdf"):
    fig.savefig(f"fig_central_illustration.{ext}", dpi=300)
print("wrote fig_central_illustration.png / .pdf")
