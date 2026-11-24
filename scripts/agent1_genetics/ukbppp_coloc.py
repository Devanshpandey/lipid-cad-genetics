#!/usr/bin/env python3
"""UKB-PPP (Olink) cis-pQTL x UKB GWAS colocalization -- cross-platform check on deCODE.

BUILD TRAP. UKB-PPP rows carry GENPOS in GRCh38 but embed the GRCh37 position in the
ID string: "1:55505647:G:T:imp:v1" has GENPOS 55039974. The project GWAS are GRCh37
with rsIDs. So we merge on the GRCh37 coordinate parsed OUT OF THE ID, never on GENPOS.
Matching GENPOS against a GRCh37 window would misalign by ~466 kb and silently return
an empty or wrong overlap.

coloc.abf implementation, priors and prior SDs are byte-identical to the deCODE run
(cispqtl_pipeline.py) so the two platforms are directly comparable.
"""
import os, io, gzip, math, glob, subprocess
import numpy as np, pandas as pd

PPP    = os.environ.get("PPPDIR", "/path/to/scratch/ukbppp")
GWASD  = "/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/merged"
OUT    = os.environ.get("OUTDIR", "/path/to/scratch/ukbppp")
P1, P2, P12 = 1e-4, 1e-4, 1e-5
W_QUANT, W_CC = 0.15, 0.20
CC = {"CAD", "MI"}

# GRCh37 cis-windows, identical to the deCODE run's h19 coordinates
GENES = {
 "PCSK9":   dict(chrom=1,  lo=55205000,  hi=55830000,  lipid="LDL_C"),
 "APOB":    dict(chrom=2,  lo=20900000,  hi=21570000,  lipid="APOB"),
 "ANGPTL3": dict(chrom=1,  lo=62760000,  hi=63370000,  lipid="TRIGLY"),
}
GWAS = {"LDL_C":"lipids_LDL_C.txt.gz","APOB":"lipids_APOB.txt.gz","TRIGLY":"lipids_TRIGLY.txt.gz",
        "CAD":"outcomes_CAD.txt.gz","MI":"outcomes_MI.txt.gz"}

def read_ppp(gene, meta):
    hits = glob.glob(os.path.join(PPP, gene, "*", f"discovery_chr{meta['chrom']}_*.gz"))
    assert hits, f"no chr{meta['chrom']} file for {gene}"
    rows = []
    with gzip.open(hits[0], "rt") as f:
        hdr = f.readline().split()
        ix = {c: i for i, c in enumerate(hdr)}
        for line in f:
            p = line.split()
            vid = p[ix["ID"]]
            bits = vid.split(":")
            if len(bits) < 4: continue
            try: pos37 = int(bits[1])
            except ValueError: continue
            if pos37 < meta["lo"] or pos37 > meta["hi"]: continue
            try:
                beta = float(p[ix["BETA"]]); se = float(p[ix["SE"]])
                lp = float(p[ix["LOG10P"]]); info = float(p[ix["INFO"]]); af = float(p[ix["A1FREQ"]])
            except ValueError: continue
            if se <= 0 or info < 0.8 or min(af, 1-af) < 0.01: continue
            a0, a1 = bits[2].upper(), bits[3].upper()
            if len(a0) != 1 or len(a1) != 1: continue
            rows.append((pos37, a0, a1, beta, se, lp))
    d = pd.DataFrame(rows, columns=["pos37","a0","a1","beta","se","mlog10p"])
    d["key"] = d.pos37.astype(str) + "_" + np.minimum(d.a0, d.a1) + "_" + np.maximum(d.a0, d.a1)
    return d.sort_values("mlog10p", ascending=False).drop_duplicates("key").reset_index(drop=True)

def read_gwas(trait, chrom, lo, hi):
    f = os.path.join(GWASD, GWAS[trait])
    prog = (f"zcat {f} | awk 'NR==1{{next}} ($1=={chrom} && $2>={lo} && $2<={hi})"
            f"{{print $2,$4,$5,$9,$10,$12}}'")
    out = subprocess.run(["bash","-c",prog], capture_output=True, text=True).stdout
    d = pd.read_csv(io.StringIO(out), sep=" ", header=None,
                    names=["pos37","a0","a1","beta","se","log10p"])
    d = d[d.se > 0].copy()
    d["a0"] = d.a0.str.upper(); d["a1"] = d.a1.str.upper()
    d = d[(d.a0.str.len()==1) & (d.a1.str.len()==1)]
    d["key"] = d.pos37.astype(str) + "_" + np.minimum(d.a0, d.a1) + "_" + np.maximum(d.a0, d.a1)
    return d.sort_values("log10p", ascending=False).drop_duplicates("key").reset_index(drop=True)

def _lse(a):
    a = np.asarray(a, float); m = a.max(); return m + math.log(np.exp(a-m).sum())
def labf(beta, se, W):
    V = se*se; r = (W*W)/(W*W+V); z2 = (beta*beta)/V
    return 0.5*(np.log(1.0-r) + r*z2)
def coloc_abf(d1, d2, W1, W2):
    m = d1.merge(d2, on="key", suffixes=("_1","_2"))
    if len(m) == 0: return dict(nsnps=0, H3=float("nan"), H4=float("nan"), lead=None)
    l1 = labf(m.beta_1.values, m.se_1.values, W1)
    l2 = labf(m.beta_2.values, m.se_2.values, W2)
    L1, L2, L4 = _lse(l1), _lse(l2), _lse(l1+l2)
    s12 = L1 + L2
    lH3 = -math.inf if L4 >= s12 else s12 + math.log1p(-math.exp(L4 - s12))
    lH = [0.0, math.log(P1)+L1, math.log(P2)+L2, math.log(P1)+math.log(P2)+lH3, math.log(P12)+L4]
    den = _lse(lH); pp = [math.exp(x-den) for x in lH]
    return dict(nsnps=len(m), H3=pp[3], H4=pp[4], lead=m.loc[(l1+l2).argmax(), "key"])

rows = []
for g, meta in GENES.items():
    dp = read_ppp(g, meta)
    print(f"[{g}] {len(dp)} Olink cis variants (GRCh37 window), "
          f"lead {dp.iloc[0].key} log10P={dp.iloc[0].mlog10p:.1f}", flush=True)
    for trait in [meta["lipid"], "CAD", "MI"]:
        dg = read_gwas(trait, meta["chrom"], meta["lo"], meta["hi"])
        r = coloc_abf(dp[["key","beta","se"]], dg[["key","beta","se"]],
                      W_QUANT, W_CC if trait in CC else W_QUANT)
        rows.append(dict(gene=g, platform="Olink_UKBPPP", protein_lead=dp.iloc[0].key,
                         protein_mlog10p=round(float(dp.iloc[0].mlog10p),1),
                         comparison=f"pQTL_vs_{trait}", nsnps=r["nsnps"],
                         PP_H3=round(r["H3"],3), PP_H4=round(r["H4"],3), coloc_lead=r["lead"]))
        print(f"   vs {trait:<7} n={r['nsnps']:<6} PP.H4={r['H4']:.3f}", flush=True)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "ukbppp_coloc_results.csv"), index=False)
print("\n" + res.to_string(index=False))
