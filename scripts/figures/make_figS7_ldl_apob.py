#!/usr/bin/env python3
"""Supplementary Figure S7 - LDL-C and ApoB are genetically inseparable.
(a) MR-BMA marginal inclusion probabilities (default prior + range across priors).
(b) ApoB vs LDL-C instrument effects: 94% shared, 19/345 discordant loci."""
import os, sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, panel, title_block, save, despine  # noqa
HERE=os.path.dirname(os.path.abspath(__file__)); SUP=os.path.join(HERE,"..","supplementary")

fig=plt.figure(figsize=(10.4,4.3)); gs=fig.add_gridspec(1,2,width_ratios=[1,1.05],wspace=0.32)
axA,axB=fig.add_subplot(gs[0]),fig.add_subplot(gs[1])

# ---- (a) MR-BMA MIP ----
m=pd.read_csv(os.path.join(SUP,"TableS17_mrbma_mip.csv"))
cols={"MIP_Lpa":"Lp(a)","MIP_TG":"Triglycerides","MIP_LDL_C":"LDL-C","MIP_ApoA1":"ApoA1","MIP_HDL_C":"HDL-C","MIP_ApoB":"ApoB"}
dflt=m[(m.prior_sigma==0.5)&(m.prior_inclusion==0.5)].iloc[0]
rows=[(lab, dflt[c], m[c].min(), m[c].max()) for c,lab in cols.items()]
rows=sorted(rows,key=lambda r:r[1])
y=np.arange(len(rows))
for i,(lab,v,lo,hi) in enumerate(rows):
    sub = lab in ("LDL-C","ApoB")
    c = CB["red"] if sub else CB["green"]
    axA.barh(i, v, color=c, alpha=0.9 if sub else 0.55, height=0.62, zorder=2)
    axA.plot([lo,hi],[i,i],color=INK,lw=1.1,zorder=3)                      # prior range
    axA.text(v+0.02 if v<0.9 else v-0.02, i, f"{v:.2f}", va="center",
             ha="left" if v<0.9 else "right", fontsize=6.8, color=INK if v<0.9 else "white")
axA.set_yticks(y); axA.set_yticklabels([r[0] for r in rows], fontsize=8)
axA.set_xlim(0,1.08); axA.set_xlabel("marginal inclusion probability (MR-BMA)")
axA.axvline(0.5,color=GREY,ls=":",lw=0.8)
axA.text(0.02,len(rows)-0.4,"LDL-C + ApoB MIP $\\approx$ 1.0\ncombined: statistical substitutes",
         fontsize=6.8,color=CB["red"],va="top",fontweight="bold")
despine(axA); panel(axA,"a")
title_block(axA,"MR-BMA cannot separate LDL-C and ApoB","bars = default prior; whiskers = range across 6 priors")

# ---- (b) discordant scatter ----
d=pd.read_csv(os.path.join(HERE,"data_discordant_apob_ldl.csv"))
disc=d.stdres.abs()>2
axB.axline((0,0),slope=1,color=GREY,ls=(0,(4,3)),lw=0.9,zorder=1)
axB.scatter(d.bLDL[~disc], d.bApoB[~disc], s=14, color=GREY, alpha=0.55, edgecolor="none", zorder=2, label="concordant")
axB.scatter(d.bLDL[disc],  d.bApoB[disc],  s=26, color=CB["red"], edgecolor="white", lw=0.4, zorder=3, label=f"discordant (n={disc.sum()})")
r=np.corrcoef(d.bLDL,d.bApoB)[0,1]
axB.set_xlabel("LDL-C effect per allele (SD)"); axB.set_ylabel("ApoB effect per allele (SD)")
axB.text(0.04,0.94,f"$r$ = {r:.2f}   $R^2$ = {r**2:.2f}\n{disc.sum()}/{len(d)} discordant loci",
         transform=axB.transAxes,fontsize=7.2,va="top",color=INK)
axB.legend(loc="lower right",fontsize=7); despine(axB); panel(axB,"b")
title_block(axB,"ApoB and LDL-C effects are 94% shared","only ~6% independent variation")
save(fig,"figS7_ldl_apob")
