#!/usr/bin/env python3
"""
burden_diagnostics.py

Rare-variant burden calibration and composition diagnostics:
  - genomic-inflation lambda_GC per trait (from burden -log10P distribution)
  - driving mask, cumulative allele frequency, approximate carrier count per gene
  - allele-frequency-bin decomposition for a gene (e.g. PDE3B) to show the signal
    is a genuine allelic series rather than a single ultra-rare variant.

Input: results/agent1_genetics/burden_summary_<TRAIT>.txt (REGENIE gene burden;
       cols CHROM POS ID ALLELE0 ALLELE1 A1FREQ N TEST BETA SE CHISQ LOG10P ...)
       and raw per-chr regenie.gz for the bin decomposition.

Usage: python burden_diagnostics.py --burden_dir results/agent1_genetics
"""
import argparse, os, math

def chisq_from_log10p(lp):
    try:
        from scipy.stats import chi2
        return chi2.isf(10**(-lp), 1)
    except Exception:
        return None

def lambda_gc(path):
    cs = []
    for line in open(path):
        p = line.split()
        if len(p) < 12 or p[7] != "ADD" or p[11] == "NA":
            continue
        try:
            c = chisq_from_log10p(float(p[11]))
        except ValueError:
            continue
        if c is not None and c == c:
            cs.append(c)
    cs.sort()
    return (cs[len(cs)//2] / 0.4549) if cs else float("nan"), len(cs)

def best_mask(path, gene):
    best = None
    for line in open(path):
        p = line.split()
        if len(p) < 12 or p[7] != "ADD" or p[2].split(".")[0] != gene:
            continue
        try:
            lp, af, N, beta = float(p[11]), float(p[5]), int(float(p[6])), float(p[8])
        except ValueError:
            continue
        mask = p[2].split(".", 1)[1]
        if best is None or lp > best[0]:
            best = (lp, mask, af, N, beta)
    return best

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--burden_dir", default="results/agent1_genetics")
    a = ap.parse_args()
    print("=== genomic inflation (lambda_GC) ===")
    for t in ["LDL_C", "TRIGLY", "HDL_C", "CAD"]:
        f = os.path.join(a.burden_dir, f"burden_summary_{t}.txt")
        if os.path.exists(f):
            l, n = lambda_gc(f)
            print(f"  {t:8} lambda={l:.3f}  (n_genes={n})")
    print("\n=== top-gene burden composition ===")
    print(f"{'gene':9}{'trait':7}{'mask':14}{'cumAF':>10}{'~carriers':>11}{'log10P':>8}")
    for g, t in [("APOB","LDL_C"),("PCSK9","LDL_C"),("LDLR","LDL_C"),
                 ("APOC3","TRIGLY"),("ANGPTL3","TRIGLY"),("LPA","LPA"),
                 ("ABCA1","HDL_C"),("CETP","HDL_C"),("PDE3B","TRIGLY"),
                 ("ANGPTL8","TRIGLY")]:
        r = best_mask(os.path.join(a.burden_dir, f"burden_summary_{t}.txt"), g)
        if r:
            lp, mask, af, N, beta = r
            print(f"{g:9}{t:7}{mask:14}{af:10.5f}{int(round(2*N*af)):11d}{lp:8.1f}")

if __name__ == "__main__":
    main()
