#!/usr/bin/env python3
"""
Agent 2 — Step 1: Positional gene mapping
Maps each credible set variant to nearest gene(s) within 500kb using
Ensembl gene coordinates downloaded from UCSC (hg38).

Inputs:
  --credible_sets   finemapping_credible_sets.csv
  --coloc_results   coloc_results.csv
  --gene_bed        Ensembl genes BED file (hg38, auto-downloaded if absent)
  --out             agent2_genes/positional_mapping.csv

Usage:
  python 01_positional_mapping.py \
    --credible_sets /path/to/finemapping_credible_sets.csv \
    --coloc_results /path/to/coloc_results.csv \
    --gene_bed      /path/to/reference/ensembl_genes_hg38.bed \
    --out           /path/to/agent2_genes/positional_mapping.csv
"""

import argparse, os, sys, subprocess
import pandas as pd
import numpy as np

random = 42

def download_gene_bed(out_path):
    """Download Ensembl gene coordinates from UCSC for hg38."""
    print(f"  Downloading Ensembl gene coordinates to {out_path}...")
    url = ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ensGene.txt.gz")
    names_url = ("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ensemblToGeneName.txt.gz")

    tmp = out_path + ".tmp"
    subprocess.run(["wget", "-q", "-O", tmp + ".gz", url], check=True)
    subprocess.run(["wget", "-q", "-O", out_path + ".names.gz", names_url], check=True)

    # Parse UCSC ensGene table: bin,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,...
    genes = pd.read_csv(tmp + ".gz", sep="\t", header=None, compression="gzip",
        names=["bin","transcript","chrom","strand","txStart","txEnd",
               "cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score",
               "name2","cdsStartStat","cdsEndStat","exonFrames"])
    names = pd.read_csv(out_path + ".names.gz", sep="\t", header=None,
                        compression="gzip", names=["transcript","gene_name"])
    genes = genes.merge(names, on="transcript", how="left")
    genes["gene_name"] = genes["gene_name"].fillna(genes["name2"])
    genes["chrom"] = genes["chrom"].str.replace("chr", "").str.strip()
    # Keep protein-coding-like entries (filter random contigs)
    genes = genes[genes["chrom"].str.match(r"^\d+$|^X$|^Y$")]
    genes["chrom"] = genes["chrom"].astype(str)

    bed = genes[["chrom","txStart","txEnd","gene_name","transcript","strand"]].drop_duplicates()
    bed.to_csv(out_path, sep="\t", index=False, header=True)
    os.remove(tmp + ".gz")
    print(f"  Gene BED written: {len(bed)} transcripts")
    return out_path


def map_variants_to_genes(variants_df, genes_df, window_kb=500):
    """
    For each variant, find all genes within window_kb.
    Returns dataframe with variant-gene pairs and distance.
    """
    window = window_kb * 1000
    results = []

    for _, var in variants_df.iterrows():
        chrom = str(var["chr"])
        pos   = int(var["locus_bp"])
        snp   = var["top_snp"]
        exp   = var["exposure"]
        out   = var["outcome"]
        pip   = var["top_pip"]
        pp4   = var["PP.H4"]

        nearby = genes_df[
            (genes_df["chrom"] == chrom) &
            (genes_df["txEnd"]   >= pos - window) &
            (genes_df["txStart"] <= pos + window)
        ].copy()

        if len(nearby) == 0:
            # No gene within window — record nearest gene regardless
            same_chr = genes_df[genes_df["chrom"] == chrom].copy()
            if len(same_chr) == 0:
                continue
            same_chr["dist"] = same_chr.apply(
                lambda r: max(0, r["txStart"] - pos, pos - r["txEnd"]), axis=1)
            nearby = same_chr.nsmallest(1, "dist")

        nearby["dist_to_variant"] = nearby.apply(
            lambda r: max(0, r["txStart"] - pos, pos - r["txEnd"]), axis=1)
        nearby["within_gene"] = (
            (nearby["txStart"] <= pos) & (nearby["txEnd"] >= pos)
        ).astype(int)
        nearby["is_nearest"] = (nearby["dist_to_variant"] == nearby["dist_to_variant"].min()).astype(int)

        for _, gene in nearby.iterrows():
            results.append({
                "exposure":        exp,
                "outcome":         out,
                "top_snp":         snp,
                "chr":             chrom,
                "locus_bp":        pos,
                "top_pip":         pip,
                "PP.H4":           pp4,
                "gene_name":       gene["gene_name"],
                "transcript":      gene["transcript"],
                "gene_start":      gene["txStart"],
                "gene_end":        gene["txEnd"],
                "strand":          gene["strand"],
                "dist_to_variant": int(gene["dist_to_variant"]),
                "within_gene":     int(gene["within_gene"]),
                "is_nearest":      int(gene["is_nearest"]),
            })

    return pd.DataFrame(results)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--credible_sets", required=True)
    parser.add_argument("--coloc_results",  required=True)
    parser.add_argument("--gene_bed",       required=True)
    parser.add_argument("--window_kb",      type=int, default=500)
    parser.add_argument("--out",            required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # ── Load gene coordinates ─────────────────────────────────────
    if not os.path.exists(args.gene_bed):
        download_gene_bed(args.gene_bed)

    print("[positional] Loading gene coordinates...")
    genes = pd.read_csv(args.gene_bed, sep="\t", dtype={"chrom": str})
    genes = genes.drop_duplicates(subset=["chrom","txStart","txEnd","gene_name"])
    print(f"  {len(genes)} gene entries loaded")

    # ── Load credible sets (top SNP per CS) ───────────────────────
    print("[positional] Loading credible sets...")
    cs = pd.read_csv(args.credible_sets)
    # Deduplicate: keep one row per unique locus x exposure (highest PIP)
    cs = cs.sort_values("top_pip", ascending=False)
    cs_top = cs.drop_duplicates(subset=["exposure", "chr", "locus_bp"])
    print(f"  {len(cs_top)} unique exposure x locus combinations")

    # ── Map variants to genes ─────────────────────────────────────
    print(f"[positional] Mapping variants to genes (window={args.window_kb}kb)...")
    results = map_variants_to_genes(cs_top, genes, window_kb=args.window_kb)

    if len(results) == 0:
        print("[positional] WARNING: No gene mappings found")
        sys.exit(0)

    # ── Add positional evidence score ─────────────────────────────
    # Score: 3 if variant is within gene body, 2 if nearest gene, 1 if within window
    results["positional_score"] = results.apply(lambda r:
        3 if r["within_gene"] else (2 if r["is_nearest"] else 1), axis=1)

    # ── Aggregate to gene level (max score across loci) ──────────
    gene_summary = results.groupby("gene_name").agg(
        n_loci        = ("locus_bp", "nunique"),
        max_pip       = ("top_pip",  "max"),
        max_pp4       = ("PP.H4",    "max"),
        positional_score = ("positional_score", "max"),
        loci          = ("locus_bp", lambda x: ";".join(str(v) for v in x.unique())),
        exposures     = ("exposure", lambda x: ";".join(x.unique())),
        outcomes      = ("outcome",  lambda x: ";".join(x.unique())),
        min_dist      = ("dist_to_variant", "min"),
    ).reset_index()

    # Save full mapping and summary
    results.to_csv(args.out.replace(".csv", "_full.csv"), index=False)
    gene_summary.to_csv(args.out, index=False)

    print(f"[positional] {len(results)} variant-gene pairs across {len(gene_summary)} unique genes")
    print(f"[positional] Written: {args.out}")

    # Print top genes
    top = gene_summary.nlargest(20, "max_pip")
    print("\nTop 20 genes by max credible set PIP:")
    print(top[["gene_name","n_loci","max_pip","max_pp4","positional_score","min_dist","exposures"]].to_string(index=False))


if __name__ == "__main__":
    main()
