#!/usr/bin/env python3
"""
Agent 2 — Step 3: OpenTargets L2G scores + druggability
Queries OpenTargets Genetics locus-to-gene (L2G) scores for our
colocalized loci, and OpenTargets Platform for drug target annotation.

Inputs:
  --credible_sets   finemapping_credible_sets.csv
  --ot_l2g_dir      /path/to/opentargets/l2g/     (parquet files)
  --ot_target_dir   /path/to/opentargets/target/   (parquet files)
  --ot_drug_dir     /path/to/opentargets/known_drug/
  --out             agent2_genes/opentargets_l2g.csv

Usage:
  python 03_opentargets_l2g.py \
    --credible_sets /path/to/finemapping_credible_sets.csv \
    --ot_l2g_dir    /path/to/reference/opentargets/l2g \
    --ot_target_dir /path/to/reference/opentargets/target \
    --ot_drug_dir   /path/to/reference/opentargets/known_drug \
    --out           /path/to/agent2_genes/opentargets_l2g.csv
"""

import argparse, os, sys, glob
import pandas as pd
import numpy as np

# Known priority genes from CLAUDE.md — pre-annotated
PRIORITY_GENES = {
    # Lipid metabolism
    "LDLR":    {"category": "lipid_metabolism",  "druggable": True},
    "APOB":    {"category": "lipid_metabolism",  "druggable": True},
    "PCSK9":   {"category": "lipid_metabolism",  "druggable": True},
    "APOE":    {"category": "lipid_metabolism",  "druggable": False},
    "LPA":     {"category": "lipid_metabolism",  "druggable": True},
    "SORT1":   {"category": "lipid_metabolism",  "druggable": False},
    "ANGPTL3": {"category": "lipid_metabolism",  "druggable": True},
    "APOC3":   {"category": "lipid_metabolism",  "druggable": True},
    "TRIB1":   {"category": "lipid_metabolism",  "druggable": True},
    "MLXIPL":  {"category": "lipid_metabolism",  "druggable": False},
    "LPL":     {"category": "lipid_metabolism",  "druggable": True},
    # Inflammation
    "IL6R":    {"category": "inflammation",      "druggable": True},
    "NLRP3":   {"category": "inflammation",      "druggable": True},
    "CXCL12":  {"category": "inflammation",      "druggable": False},
    # Plaque biology
    "PHACTR1": {"category": "plaque_biology",    "druggable": False},
    "ADAMTS7": {"category": "plaque_biology",    "druggable": True},
    "COL4A1":  {"category": "plaque_biology",    "druggable": False},
    "COL4A2":  {"category": "plaque_biology",    "druggable": False},
    "TCF21":   {"category": "plaque_biology",    "druggable": False},
    # Thrombosis
    "F5":      {"category": "thrombosis",        "druggable": True},
    "F2":      {"category": "thrombosis",        "druggable": True},
    "ABO":     {"category": "thrombosis",        "druggable": False},
    "VWF":     {"category": "thrombosis",        "druggable": True},
    # BP / remodeling
    "NOS3":    {"category": "bp_remodeling",     "druggable": True},
    "NPR3":    {"category": "bp_remodeling",     "druggable": True},
    "CYP17A1": {"category": "bp_remodeling",     "druggable": True},
}


def load_parquet_dir(directory, columns=None):
    """Load all parquet files in a directory into a single dataframe."""
    files = glob.glob(os.path.join(directory, "*.parquet"))
    if not files:
        files = glob.glob(os.path.join(directory, "**/*.parquet"), recursive=True)
    if not files:
        print(f"  WARNING: No parquet files in {directory}")
        return pd.DataFrame()
    dfs = []
    for f in files:
        try:
            df = pd.read_parquet(f, columns=columns) if columns else pd.read_parquet(f)
            dfs.append(df)
        except Exception as e:
            print(f"  Warning: could not read {f}: {e}")
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


def get_l2g_scores(cs_df, l2g_dir, window_mb=1.0):
    """
    Query L2G parquet files for genes near our loci.
    L2G schema: study_id, variant_id, gene_id, gene_name, y_proba_full_model, ...
    We filter by chromosome position proximity.
    """
    if not os.path.exists(l2g_dir):
        print(f"  L2G directory not found: {l2g_dir}")
        return pd.DataFrame()

    print(f"  Loading OpenTargets L2G scores...")
    l2g = load_parquet_dir(l2g_dir, columns=None)
    if l2g.empty:
        return pd.DataFrame()

    print(f"  L2G table: {len(l2g)} rows, columns: {list(l2g.columns[:10])}")

    # L2G has variant positions encoded in variant_id (e.g. "1_55505647_C_T")
    # Extract chrom and pos
    if "variant_id" in l2g.columns:
        parts = l2g["variant_id"].str.split("_", expand=True)
        l2g["chrom_l2g"] = parts[0].astype(str)
        l2g["pos_l2g"]   = pd.to_numeric(parts[1], errors="coerce")
    elif "chrom" in l2g.columns and "pos" in l2g.columns:
        l2g["chrom_l2g"] = l2g["chrom"].astype(str)
        l2g["pos_l2g"]   = l2g["pos"]

    results = []
    for _, row in cs_df.iterrows():
        chrom = str(row["chr"])
        bp    = int(row["locus_bp"])
        exp   = row["exposure"]
        snp   = row["top_snp"]

        window = int(window_mb * 1e6)
        nearby = l2g[
            (l2g["chrom_l2g"] == chrom) &
            (l2g["pos_l2g"] >= bp - window) &
            (l2g["pos_l2g"] <= bp + window)
        ].copy()

        if nearby.empty:
            continue

        # Score column name varies by OT version
        score_col = next((c for c in ["y_proba_full_model", "score", "l2g_score"]
                          if c in nearby.columns), None)
        gene_col  = next((c for c in ["gene_name", "gene_symbol", "target_id"]
                          if c in nearby.columns), None)

        if score_col is None or gene_col is None:
            continue

        nearby = nearby.sort_values(score_col, ascending=False)

        for _, g in nearby.iterrows():
            results.append({
                "exposure":    exp,
                "locus_chr":   chrom,
                "locus_bp":    bp,
                "top_snp":     snp,
                "gene_name":   g.get(gene_col, ""),
                "l2g_score":   g.get(score_col, np.nan),
                "PP.H4_gwas":  row.get("PP.H4", np.nan),
            })

    return pd.DataFrame(results)


def get_drug_targets(gene_list, target_dir, drug_dir):
    """
    Flag genes as drug targets using OpenTargets Platform.
    Returns dict: gene_name -> {"druggable": bool, "approved_drugs": list, "tractability": str}
    """
    result = {}

    # Load target info
    if os.path.exists(target_dir):
        targets = load_parquet_dir(target_dir, columns=None)
        if not targets.empty:
            name_col = next((c for c in ["approvedSymbol","gene_name","symbol"]
                             if c in targets.columns), None)
            tract_col = next((c for c in ["tractability","hasTractabilityData"]
                              if c in targets.columns), None)
            if name_col:
                for _, row in targets[targets[name_col].isin(gene_list)].iterrows():
                    gname = row[name_col]
                    result[gname] = {
                        "druggable": True,
                        "tractability": str(row.get(tract_col, "")) if tract_col else "",
                        "approved_drugs": []
                    }

    # Load known drugs
    if os.path.exists(drug_dir):
        drugs = load_parquet_dir(drug_dir, columns=None)
        if not drugs.empty:
            gene_col = next((c for c in ["targetSymbol","gene_name","symbol"]
                             if c in drugs.columns), None)
            drug_col = next((c for c in ["drugId","drug_name","label"]
                             if c in drugs.columns), None)
            if gene_col and drug_col:
                for gene, grp in drugs[drugs[gene_col].isin(gene_list)].groupby(gene_col):
                    if gene not in result:
                        result[gene] = {"druggable": True, "tractability": "", "approved_drugs": []}
                    result[gene]["approved_drugs"] = list(grp[drug_col].unique())

    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--credible_sets",  required=True)
    parser.add_argument("--ot_l2g_dir",     required=True)
    parser.add_argument("--ot_target_dir",  required=True)
    parser.add_argument("--ot_drug_dir",    required=True)
    parser.add_argument("--out",            required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # ── Load credible sets ─────────────────────────────────────────
    print("[OpenTargets] Loading credible sets...")
    cs = pd.read_csv(args.credible_sets)
    cs = cs.sort_values("top_pip", ascending=False).drop_duplicates(
        subset=["exposure", "chr", "locus_bp"])
    print(f"  {len(cs)} unique loci")

    # ── L2G scores ────────────────────────────────────────────────
    l2g_df = get_l2g_scores(cs, args.ot_l2g_dir)

    # ── Drug target annotation ────────────────────────────────────
    all_genes = list(l2g_df["gene_name"].unique()) if not l2g_df.empty else []
    all_genes += list(PRIORITY_GENES.keys())
    drug_info = get_drug_targets(list(set(all_genes)), args.ot_target_dir, args.ot_drug_dir)

    # ── Merge with priority gene annotation ───────────────────────
    if not l2g_df.empty:
        l2g_df["l2g_sig"]         = l2g_df["l2g_score"] >= 0.5
        l2g_df["is_priority_gene"] = l2g_df["gene_name"].isin(PRIORITY_GENES)
        l2g_df["gene_category"]    = l2g_df["gene_name"].map(
            lambda g: PRIORITY_GENES.get(g, {}).get("category", "other"))
        l2g_df["druggable_ot"]     = l2g_df["gene_name"].map(
            lambda g: drug_info.get(g, {}).get("druggable",
                PRIORITY_GENES.get(g, {}).get("druggable", False)))
        l2g_df["approved_drugs"]   = l2g_df["gene_name"].map(
            lambda g: ";".join(drug_info.get(g, {}).get("approved_drugs", [])))

        # L2G evidence score: 3 if >= 0.8, 2 if >= 0.5, 1 if >= 0.2
        l2g_df["l2g_score_pts"] = pd.cut(l2g_df["l2g_score"],
            bins=[-1, 0.2, 0.5, 0.8, 2],
            labels=[0, 1, 2, 3]).astype(float)

        l2g_df.to_csv(args.out, index=False)
        n_sig = l2g_df["l2g_sig"].sum()
        print(f"\n[OpenTargets] {len(l2g_df)} gene-locus pairs, {n_sig} with L2G >= 0.5")
        print(f"[OpenTargets] Written: {args.out}")

        top = l2g_df[l2g_df["l2g_sig"]].nlargest(15, "l2g_score")
        if len(top) > 0:
            print("\nTop L2G genes:")
            print(top[["gene_name","exposure","locus_chr","locus_bp",
                       "l2g_score","gene_category","druggable_ot"]].to_string(index=False))
    else:
        print("[OpenTargets] WARNING: No L2G data found")
        # Still write priority gene annotations
        priority_df = pd.DataFrame([
            {"gene_name": g, **info, "l2g_score": np.nan, "l2g_sig": False,
             "approved_drugs": ";".join(drug_info.get(g, {}).get("approved_drugs", []))}
            for g, info in PRIORITY_GENES.items()
        ])
        priority_df.to_csv(args.out, index=False)
        print(f"[OpenTargets] Written priority gene annotations: {args.out}")


if __name__ == "__main__":
    main()
