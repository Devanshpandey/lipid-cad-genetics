#!/usr/bin/env python3
"""deCODE (SomaScan v4) cis-pQTL x UKB GWAS colocalization -- SomaScan arm.

BUILD. deCODE Pos is GRCh38; the project GWAS are GRCh37. We never compare the
two coordinate systems. Instead we take the rsIDs inside the GRCh37 cis window
from the GWAS and select deCODE rows by rsID membership, which is build-agnostic.
(The Olink arm has the mirror-image trap: its GENPOS is GRCh38 while its ID
string carries GRCh37, so there we parse the position out of the ID.)

coloc.abf implementation, priors and prior SDs are identical to the Olink run
(ukbppp_coloc_v2.py) so the two platforms are directly comparable -- which is
the entire point, given Eldjarn 2023 reports cis-pQTL assay support for 43% of
SomaScan versus 72% of Olink assays.
"""
import os, io, gzip, math, subprocess
import numpy as np, pandas as pd

DEC   = os.environ.get("DECDIR", "/path/to/scratch/pqtl2/decode")
DEC2  = os.environ.get("DECDIR2", "/path/to/scratch/cispqtl")
GWASD = "/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/merged"
OUT   = os.environ.get("OUTDIR", "/path/to/scratch/pqtl2")
P1, P2, P12 = 1e-4, 1e-4, 1e-5
W_QUANT, W_CC = 0.15, 0.20
CC = {"CAD", "MI"}

# GRCh37 cis-windows, identical to the Olink arm so the platforms are comparable
GENES = {
 "PCSK9":   dict(chrom=1,  lo=55205000,  hi=55830000,  lipid="LDL_C"),
 "APOB":    dict(chrom=2,  lo=20900000,  hi=21570000,  lipid="APOB"),
 "ANGPTL3": dict(chrom=1,  lo=62760000,  hi=63370000,  lipid="TRIGLY"),
 "APOC3":   dict(chrom=11, lo=116400000, hi=117010000, lipid="TRIGLY"),
 "LDLR":    dict(chrom=19, lo=10900000,  hi=11545000,  lipid="LDL_C"),
 "SORT1":   dict(chrom=1,  lo=109550000, hi=110240000, lipid="LDL_C"),
 "APOE":    dict(chrom=19, lo=45100000,  hi=45715000,  lipid="LDL_C"),
 "APOC1":   dict(chrom=19, lo=45115000,  hi=45725000,  lipid="LDL_C"),
 "APOA5":   dict(chrom=11, lo=116490000, hi=117094000, lipid="TRIGLY"),
 "HMGCR":   dict(chrom=5,  lo=74333000,  hi=74958000,  lipid="LDL_C"),
 "LDLRAP1": dict(chrom=1,  lo=25570000,  hi=26195000,  lipid="LDL_C"),
 "APOA1":   dict(chrom=11, lo=116405000, hi=117010000, lipid="APOA1"),
 "ANGPTL4": dict(chrom=19, lo=8125000,   hi=8740000,   lipid="TRIGLY"),
 "LCAT":    dict(chrom=16, lo=67670000,  hi=68280000,  lipid="HDL_C"),
}
GWAS = {"LDL_C":"lipids_LDL_C.txt.gz","APOB":"lipids_APOB.txt.gz","TRIGLY":"lipids_TRIGLY.txt.gz",
        "HDL_C":"lipids_HDL_C.txt.gz","APOA1":"lipids_APOA1.txt.gz",
        "CAD":"outcomes_CAD.txt.gz","MI":"outcomes_MI.txt.gz"}

def gwas_window(trait, chrom, lo, hi):
    f = os.path.join(GWASD, GWAS[trait])
    prog = (f"zcat {f} | awk 'NR==1{{next}} ($1=={chrom} && $2>={lo} && $2<={hi})"
            f"{{print $3,$9,$10,$12}}'")
    out = subprocess.run(["bash","-c",prog], capture_output=True, text=True).stdout
    d = pd.read_csv(io.StringIO(out), sep=" ", header=None, names=["rsid","beta","se","log10p"])
    d = d[(d.se > 0) & d.rsid.astype(str).str.startswith("rs")]
    return d.sort_values("log10p", ascending=False).drop_duplicates("rsid").reset_index(drop=True)

def decode_window(gene, chrom, rsids):
    path = None
    for cand in (os.path.join(DEC, f"{gene}.full.txt.gz"), os.path.join(DEC2, f"{gene}.full.txt.gz")):
        if os.path.exists(cand): path = cand; break
    if path is None: return None
    keep = set(rsids)
    rows = []
    want = f"chr{chrom}".encode()
    with gzip.open(path, "rb") as f:
        f.readline()
        seen = False
        for line in f:
            p = line.split(b"\t")
            if p[0] != want:
                if seen: break          # left the target chromosome block
                continue
            seen = True
            rs = p[3].decode()
            if rs not in keep: continue
            try:
                beta = float(p[6]); se = float(p[9]); lp = float(p[8]); maf = float(p[11])
            except (ValueError, IndexError): continue
            if se <= 0 or maf < 0.01: continue
            rows.append((rs, beta, se, lp))
    d = pd.DataFrame(rows, columns=["rsid","beta","se","mlog10p"])
    if len(d) == 0: return d
    return d.sort_values("mlog10p", ascending=False).drop_duplicates("rsid").reset_index(drop=True)

def _lse(a):
    a = np.asarray(a, float); m = a.max(); return m + math.log(np.exp(a-m).sum())
def labf(beta, se, W):
    V = se*se; r = (W*W)/(W*W+V); z2 = (beta*beta)/V
    return 0.5*(np.log(1.0-r) + r*z2)
def coloc_abf(d1, d2, W1, W2):
    m = d1.merge(d2, on="rsid", suffixes=("_1","_2"))
    if len(m) == 0: return dict(nsnps=0, H3=float("nan"), H4=float("nan"), lead=None)
    l1 = labf(m.beta_1.values, m.se_1.values, W1)
    l2 = labf(m.beta_2.values, m.se_2.values, W2)
    L1, L2, L4 = _lse(l1), _lse(l2), _lse(l1+l2)
    s12 = L1 + L2
    lH3 = -math.inf if L4 >= s12 else s12 + math.log1p(-math.exp(L4 - s12))
    lH = [0.0, math.log(P1)+L1, math.log(P2)+L2, math.log(P1)+math.log(P2)+lH3, math.log(P12)+L4]
    den = _lse(lH); pp = [math.exp(x-den) for x in lH]
    return dict(nsnps=len(m), H3=pp[3], H4=pp[4], lead=m.loc[(l1+l2).argmax(), "rsid"])

rows = []
for g, meta in GENES.items():
    anchor = gwas_window(meta["lipid"], meta["chrom"], meta["lo"], meta["hi"])
    dp = decode_window(g, meta["chrom"], anchor.rsid.tolist())
    if dp is None:
        print(f"[{g}] no deCODE file - skipping", flush=True); continue
    if len(dp) == 0:
        print(f"[{g}] no cis variants after rsID match - skipping", flush=True); continue
    print(f"[{g}] {len(dp)} SomaScan cis variants, lead {dp.iloc[0].rsid} "
          f"log10P={dp.iloc[0].mlog10p:.1f}", flush=True)
    for trait in [meta["lipid"], "CAD", "MI"]:
        dg = anchor if trait == meta["lipid"] else gwas_window(trait, meta["chrom"], meta["lo"], meta["hi"])
        r = coloc_abf(dp[["rsid","beta","se"]], dg[["rsid","beta","se"]],
                      W_QUANT, W_CC if trait in CC else W_QUANT)
        rows.append(dict(gene=g, platform="SomaScan_deCODE", protein_lead=dp.iloc[0].rsid,
                         protein_mlog10p=round(float(dp.iloc[0].mlog10p),1),
                         comparison=f"pQTL_vs_{trait}", nsnps=r["nsnps"],
                         PP_H3=round(r["H3"],3), PP_H4=round(r["H4"],3), coloc_lead=r["lead"]))
        print(f"   vs {trait:<7} n={r['nsnps']:<6} PP.H4={r['H4']:.3f}", flush=True)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "decode_coloc_results_v2.csv"), index=False)
print("\n" + res.to_string(index=False))
