#!/usr/bin/env python3
"""Supplementary Figure S8 - Robustness and complementarity of the convergence score.
(a) Orthogonal to Open Targets L2G. (b) Layer ablation. (c) Weight sensitivity (PDE3B rank).
(d) Recovery by target class."""
import os, sys, numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import CB, INK, GREY, panel, title_block, save, despine  # noqa
HERE=os.path.dirname(os.path.abspath(__file__)); SUP=os.path.join(HERE,"..","supplementary")

fig=plt.figure(figsize=(10.6,8.0)); gs=fig.add_gridspec(2,2,hspace=0.55,wspace=0.34)
axA,axB,axC,axD=[fig.add_subplot(gs[i]) for i in range(4)]

# ---- (a) L2G orthogonality ----
d=pd.read_csv(os.path.join(SUP,"TableS19_l2g_benchmark.csv"))
rho=spearmanr(d.integrated,d.l2g).correlation
can=d.canonical==1
axA.scatter(d.l2g[~can], d.integrated[~can], s=16, color=GREY, alpha=0.55, edgecolor="none", label="other genes")
axA.scatter(d.l2g[can],  d.integrated[can],  s=42, color=CB["amber"], edgecolor="black", lw=0.4, label="canonical target")
axA.set_xlabel("Open Targets L2G score"); axA.set_ylabel("convergence (integrated) score")
axA.text(0.04,0.95,f"Spearman $\\rho$ = {rho:.2f}\n22/27 canonical targets\noutside L2G-scored set",
         transform=axA.transAxes,fontsize=7.2,va="top",color=INK)
axA.legend(loc="lower right",fontsize=6.8); despine(axA); panel(axA,"a")
title_block(axA,"Orthogonal to common-variant L2G","captures rare-variant genes L2G cannot score")

# ---- (b) layer ablation ----
ab=pd.read_csv(os.path.join(SUP,"TableS20_layer_ablation.csv"))
labs=["burden\nonly","+ coloc","+ outcome\nburden","full"]; auc=ab.AUROC_canonical_recovery.values
cols=[GREY,CB["green"],GREY,CB["green"]]
x=np.arange(4); axB.bar(x,auc,color=cols,width=0.64,zorder=2)
for xi,v in zip(x,auc): axB.text(xi,v+0.006,f"{v:.2f}",ha="center",fontsize=7,color=INK)
axB.set_xticks(x); axB.set_xticklabels(labs,fontsize=7.4); axB.set_ylim(0.8,1.0)
axB.set_ylabel("AUROC (canonical recovery)")
axB.text(0.03,0.97,"gain driven by\ncolocalization\n(+outcome: LRT $P$=0.20)",
         transform=axB.transAxes,ha="left",va="top",fontsize=6.4,color=CB["green"],fontweight="bold")
despine(axB); panel(axB,"b")
title_block(axB,"Which layer earns the gain","nested ablation of the score")

# ---- (c) weight sensitivity: PDE3B rank across specifications ----
ws=pd.read_csv(os.path.join(SUP,"TableS22_score_weight_sensitivity.csv"))
axC.hist(ws.pde3b_rank, bins=np.arange(0,26,2), color=CB["purple"], alpha=0.85, edgecolor="white")
axC.axvline(180,color=CB["red"],ls=(0,(4,3)),lw=1)
axC.text(182,axC.get_ylim()[1]*0.5,"top 1% (rank 180)",fontsize=6.6,color=CB["red"],rotation=90,va="center")
axC.set_xlabel("PDE3B genome-wide rank"); axC.set_ylabel("number of weight specifications")
axC.text(0.96,0.95,f"top 0.1% (rank 3-22)\nin {len(ws)}/{len(ws)} specifications",
         transform=axC.transAxes,ha="right",va="top",fontsize=7.2,color=INK)
despine(axC); panel(axC,"c")
title_block(axC,"Nomination robust to weighting","108 burden:coloc:outcome / cap specifications")

# ---- (d) stratified AUROC ----
st=pd.read_csv(os.path.join(SUP,"TableS23_stratified_auroc.csv"))
labs=[s.replace("approved/trial drug target","approved/trial\ndrug targets").replace("Mendelian dyslipidemia","Mendelian\ndyslipidemia").replace("all canonical","all\ncanonical") for s in st.stratum]
x=np.arange(len(st)); axD.bar(x,st.AUROC,color=[CB["amber"],CB["blue"],CB["green"]],width=0.6,zorder=2)
for xi,v,n in zip(x,st.AUROC,st.n): axD.text(xi,v+0.006,f"{v:.2f}\n(n={n})",ha="center",fontsize=6.8,color=INK)
axD.set_xticks(x); axD.set_xticklabels(labs,fontsize=7.2); axD.set_ylim(0.8,1.02)
axD.set_ylabel("AUROC"); despine(axD); panel(axD,"d")
title_block(axD,"Recovery holds across target classes","positive set analysed by evidence type")
save(fig,"figS8_score_robustness")
