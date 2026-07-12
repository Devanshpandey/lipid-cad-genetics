#!/usr/bin/env python3
"""
Supplementary Figure S12 — PDE3B carrier survival diagnostics.
  (a) Adjusted cumulative incidence of CAD in PDE3B pLoF carriers vs non-carriers
      (nearly identical), with number-at-risk.
  (b) Power to detect a carrier CAD effect given 92 events: 80% power only for
      HR <= 0.75 or >= 1.34, so the null Cox HR (1.03) is uninformative about a
      modest effect. Proportional-hazards (Schoenfeld) held (carrier P=0.10).

Outputs figS12_pde3b_survival.{pdf,png}. Data: pde3b_km.csv, pde3b_power.csv.
"""
import os, sys, csv
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, SUBTLE, LIPID, save, despine  # noqa

HERE = os.path.dirname(os.path.abspath(__file__)); DATA = os.path.join(HERE, "data")
TEAL, BURG = LIPID["rare"], LIPID["outcome"]
km = list(csv.DictReader(open(os.path.join(DATA, "pde3b_km.csv"))))
pw = list(csv.DictReader(open(os.path.join(DATA, "pde3b_power.csv"))))

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 4.7), gridspec_kw={"wspace": 0.30})

# ---------------------------------------------------- (a) cumulative incidence
for g, col, lab in (("0", GREY, "non-carrier"), ("1", BURG, "pLoF carrier")):
    t = [float(x["time"]) for x in km if x["carrier"] == g]
    c = [float(x["cuminc"]) * 100 for x in km if x["carrier"] == g]
    axA.step(t, c, where="post", color=col, lw=2.0, label=lab, zorder=3)
axA.set_xlabel("years from baseline"); axA.set_ylabel("cumulative incidence of CAD (%)")
axA.set_xlim(0, 14); axA.set_ylim(0, 10)
axA.legend(loc="upper left", fontsize=7.4); despine(axA)
# number at risk (carriers) at a few times
nr = {float(x["time"]): int(x["n_risk"]) for x in km if x["carrier"] == "1"}
ticks = [0, 4, 8, 12]
axA.text(0.0, -0.20, "carriers at risk:  " + "    ".join(f"{int(t)}y: {nr.get(float(t), '-')}" for t in ticks),
         transform=axA.transAxes, fontsize=6.0, color=SUBTLE)
axA.text(-0.02, 1.13, "a", transform=axA.transAxes, fontsize=13, fontweight="bold", va="top", ha="right", color=INK)
axA.text(0.5, 1.125, "Carrier vs non-carrier CAD incidence", transform=axA.transAxes,
         fontsize=9.4, fontweight="bold", va="top", ha="center", color=INK)
axA.text(0.5, 1.045, "curves nearly identical (Cox HR 1.03, 92 carrier events)",
         transform=axA.transAxes, fontsize=7.2, va="top", ha="center", color=SUBTLE)

# ---------------------------------------------------- (b) power curve
hr = np.array([float(x["HR"]) for x in pw]); power = np.array([float(x["power"]) for x in pw])
axB.plot(hr, power, color=CB["purple"], lw=2.2, zorder=3)
axB.axhline(0.8, color=GREY, ls=(0, (3, 3)), lw=0.9)
axB.text(1.55, 0.82, "80% power", fontsize=6.4, color=GREY, ha="right")
axB.axvspan(0.745, 1.342, color="#F3EDF7", zorder=0)
axB.text(1.043, 0.12, "under-powered\nzone (HR 0.75-1.34)", fontsize=6.2, color=CB["purple"], ha="center")
axB.axvline(1.028, color=BURG, lw=1.4, zorder=2)
axB.text(1.028, 0.55, "observed\nHR 1.03", fontsize=6.4, color=BURG, ha="left", rotation=90, va="center")
axB.set_xlabel("true carrier hazard ratio for CAD"); axB.set_ylabel("power (92 events)")
axB.set_xlim(0.5, 1.6); axB.set_ylim(0, 1.02); despine(axB)
axB.text(0.03, 0.05, "proportional hazards held\n(Schoenfeld carrier P = 0.10)",
         transform=axB.transAxes, fontsize=6.6, color=INK)
axB.text(-0.02, 1.13, "b", transform=axB.transAxes, fontsize=13, fontweight="bold", va="top", ha="right", color=INK)
axB.text(0.5, 1.125, "Underpowered to exclude a modest effect", transform=axB.transAxes,
         fontsize=9.4, fontweight="bold", va="top", ha="center", color=INK)
axB.text(0.5, 1.045, "80% power only for HR $\\leq$0.75 or $\\geq$1.34",
         transform=axB.transAxes, fontsize=7.2, va="top", ha="center", color=SUBTLE)
save(fig, "figS12_pde3b_survival")
