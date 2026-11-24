#!/usr/bin/env python3
"""Distance-pruned coloc sentinels: no LD panel required.

The LD-clumped definition (make_coloc_sentinels.R, queued as 07a/07b) needs PLINK
against the 1.4 TB genotype panel. Greedy distance pruning -- take the strongest
GWS variant, drop everything within +/-500 kb, repeat -- needs only the summary
statistics, is deterministic, and matches the +/-500 kb coloc window exactly.
It is the sensitivity analysis we would report anyway, so running it first costs
nothing and unblocks the genome-wide result while the cluster is saturated.

Thresholds stay on the LOG10P scale throughout: 10^(-LOG10P) underflows to 0 for
-log10P > 323.3, which is the bug that collapsed the original sentinel set.
"""
import gzip, sys, os
G = os.path.expandvars("/path/to/cad-genetics/results/agent1_genetics/gwas_summary_stats/merged")
THRESH = 7.30103
SEP = 500_000

for trait in ["LDL_C","HDL_C","TRIGLY","TOT_CHOL","APOA1","APOB"]:
    sig = []
    with gzip.open(f"{G}/lipids_{trait}.txt.gz", "rt") as f:
        f.readline()
        for line in f:
            p = line.split()
            try:
                lp = float(p[11])
            except (ValueError, IndexError):
                continue
            if lp > THRESH:
                sig.append((lp, int(p[0]), int(p[1]), p[2]))
    sig.sort(reverse=True)
    taken = []
    for lp, c, b, rs in sig:
        if any(tc == c and abs(tb - b) <= SEP for tc, tb in taken):
            continue
        taken.append((c, b))
        sys.stdout.write("")
    keep = []
    taken = []
    for lp, c, b, rs in sig:
        if any(tc == c and abs(tb - b) <= SEP for tc, tb in taken): continue
        taken.append((c, b)); keep.append((c, b, rs, lp))
    keep.sort()
    with open(f"sentinels_{trait}.tsv", "w") as o:
        o.write("SNP\tCHR\tBP\tLOG10P\n")
        for c, b, rs, lp in keep:
            o.write(f"{rs}\t{c}\t{b}\t{lp}\n")
    print(f"{trait}: {len(sig)} GWS variants -> {len(keep)} independent loci "
          f"on {len(set(c for c,_,_,_ in keep))} chromosomes", flush=True)
