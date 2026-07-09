#!/usr/bin/env python3
"""
integrated_convergence_score.py

Combine common- and rare-variant evidence into a single per-gene score that
credentials established lipid/CAD genes and nominates non-canonical candidates.

Inputs (produced upstream by the Agent-1 pipeline):
  results/agent1_genetics/burden_summary_<TRAIT>.txt   # REGENIE gene burden, per trait
Colocalized loci are read from the coloc step (here hard-coded to the five loci
that exceeded PP.H4 > 0.8; edit `COLOC` for a different run).

Output:
  convergence_scores.csv   # gene-ranked integrated evidence table

Usage:
  python integrated_convergence_score.py --burden_dir results/agent1_genetics --out convergence_scores.csv
"""
import argparse, os

LIPIDS   = "LDL_C HDL_C TRIGLY TOT_CHOL APOA1 APOB LPA nonHDL_C".split()
OUTCOMES = "CAD MI REVASC MACE STROKE HF CV_DEATH AF PAD".split()
COLOC    = {"PCSK9","CELSR2","SORT1","PSRC1","LDLR","APOE","APOC1","TRIB1"}
CANONICAL = {"LDLR","APOB","PCSK9","APOC3","ANGPTL3","APOA5","ABCA1","CETP","LCAT",
             "LPL","LPA","APOE","APOA1","SORT1","CELSR2","PSRC1","MTTP","NPC1L1",
             "HMGCR","PLG","ANGPTL4","APOC2","LIPC","LIPG","SCARB1","ABCG5","ABCG8"}
THR = 5.6  # exome-wide significance, -log10P (p < 2.5e-6)


def best_per_gene(path):
    """Max ADD-burden -log10P (and its beta) per gene from one trait summary."""
    best = {}
    if not os.path.exists(path):
        return best
    for line in open(path):
        p = line.split()
        if len(p) < 12 or p[7] != "ADD" or p[11] == "NA":
            continue
        try:
            lp, beta = float(p[11]), float(p[8])
        except ValueError:
            continue
        g = p[2].split(".")[0]
        if g not in best or lp > best[g][0]:
            best[g] = (lp, beta, p[0])
    return best


def aggregate(traits, burden_dir):
    agg = {}
    for t in traits:
        for g, (lp, beta, chrom) in best_per_gene(
                os.path.join(burden_dir, f"burden_summary_{t}.txt")).items():
            if g not in agg or lp > agg[g][0]:
                agg[g] = (lp, t, beta, chrom)
    return agg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--burden_dir", default="results/agent1_genetics")
    ap.add_argument("--out", default="convergence_scores.csv")
    a = ap.parse_args()

    lipbest = aggregate(LIPIDS, a.burden_dir)
    outbest = aggregate(OUTCOMES, a.burden_dir)

    rows = []
    for g, (lp, t, beta, chrom) in lipbest.items():
        if lp < THR:
            continue
        o = outbest.get(g, (0, "", 0, ""))
        coloc = 1 if g in COLOC else 0
        # integrated score: capped rare-lipid evidence + coloc bonus + pleiotropy bonus
        score = min(lp, 60) / 60 * 6 + (2 if coloc else 0) + (2 if o[0] > THR else 0)
        novel = 0 if g in CANONICAL else 1
        rows.append((round(score, 2), g, chrom, round(lp, 1), t, round(beta, 3),
                     round(o[0], 1), o[1], coloc, novel))
    rows.sort(key=lambda r: -r[0])

    with open(a.out, "w") as fh:
        fh.write("integrated_score,gene,chr,lipid_burden_log10P,lipid_trait,"
                 "lipid_beta,outcome_burden_log10P,outcome_trait,coloc_locus,"
                 "novel_noncanonical\n")
        for r in rows:
            fh.write(",".join(map(str, r)) + "\n")
    print(f"wrote {a.out}: {len(rows)} genes with exome-wide lipid burden")


if __name__ == "__main__":
    main()
