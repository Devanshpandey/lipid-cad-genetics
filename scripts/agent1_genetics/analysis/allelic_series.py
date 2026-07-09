#!/usr/bin/env python3
"""
allelic_series.py

Extract the rare-variant loss-of-function burden effect (beta, SE) on each
drug-target gene's primary lipid, for the "rare-variant allelic series
recapitulates lipid pharmacology" analysis. Pairs each gene with its matching
approved/trial therapeutic to test directional concordance.

Input:
  results/agent1_genetics/burden_summary_<TRAIT>.txt   # REGENIE gene burden
Output:
  allelic_series.csv

Usage:
  python allelic_series.py --burden_dir results/agent1_genetics --out allelic_series.csv
"""
import argparse, os

# gene -> (primary lipid trait, matching therapeutic, direction matches drug?)
TARGETS = {
    "PCSK9":   ("LDL_C",  "Evolocumab (PCSK9i) lowers LDL",        True),
    "LDLR":    ("LDL_C",  "Statins upregulate LDLR, lower LDL",    True),
    "APOB":    ("LDL_C",  "Mipomersen (APOB ASO) lowers LDL",      True),
    "APOC3":   ("TRIGLY", "Volanesorsen (APOC3 ASO) lowers TG",    True),
    "ANGPTL3": ("TRIGLY", "Evinacumab (ANGPTL3 mAb) lowers TG/LDL",True),
    "LPL":     ("TRIGLY", "LPL pathway (APOC3/ANGPTL3 inhibition)", False),
    "LPA":     ("LPA",    "Pelacarsen (LPA ASO) lowers Lp(a)",     True),
    "ABCA1":   ("HDL_C",  "Tangier gene (no approved drug)",       False),
    "CETP":    ("HDL_C",  "CETP inhibitors raise HDL",             True),
    "LCAT":    ("HDL_C",  "LCAT deficiency (no approved drug)",    False),
}


def best_burden(path, gene):
    """Strongest ADD-burden row (beta, SE, -log10P, mask) for a gene."""
    best = None
    if not os.path.exists(path):
        return best
    for line in open(path):
        p = line.split()
        if len(p) < 12 or p[7] != "ADD" or p[2].split(".")[0] != gene:
            continue
        try:
            lp, beta, se = float(p[11]), float(p[8]), float(p[9])
        except ValueError:
            continue
        mask = p[2].split(".", 1)[1]
        if best is None or lp > best[0]:
            best = (lp, beta, se, mask)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--burden_dir", default="results/agent1_genetics")
    ap.add_argument("--out", default="allelic_series.csv")
    a = ap.parse_args()

    with open(a.out, "w") as fh:
        fh.write("gene,lipid,beta_SD,se,log10P,mask,therapeutic,direction_matches_drug\n")
        for g, (trait, drug, match) in TARGETS.items():
            b = best_burden(os.path.join(a.burden_dir, f"burden_summary_{trait}.txt"), g)
            if b:
                fh.write(f"{g},{trait},{b[1]:.3f},{b[2]:.4f},{b[0]:.1f},{b[3]},"
                         f'"{drug}",{int(match)}\n')
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
