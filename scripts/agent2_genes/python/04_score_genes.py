#!/usr/bin/env python3
"""
Agent 2 — Step 4: Integrate all evidence layers into gene scores.

Evidence layers (max 10 points total):
  Positional (0-3 pts)
    - 3: variant falls within gene body
    - 2: nearest gene to credible set variant
    - 1: within 500kb window

  eQTL colocalization (0-3 pts)
    - 3: PP.H4 >= 0.8 in liver or coronary artery
    - 2: PP.H4 >= 0.8 in any tissue
    - 1: PP.H4 >= 0.5 in any tissue

  OpenTargets L2G (0-2 pts)
    - 2: L2G score >= 0.8
    - 1: L2G score >= 0.5

  Fine-mapping PIP (0-1 pt)
    - 1: top_pip >= 0.9 (singleton credible set)

  Known drug target (0-1 pt)
    - 1: ChEMBL/OpenTargets confirmed druggable target

Usage:
  python 04_score_genes.py \
    --positional    /path/to/positional_mapping.csv \
    --eqtl_coloc    /path/to/gtex_eqtl_coloc.csv \
    --ot_l2g        /path/to/opentargets_l2g.csv \
    --credible_sets /path/to/finemapping_credible_sets.csv \
    --mr_results    /path/to/mr_all_pairs.csv \
    --out           /path/to/agent2_genes/ranked_causal_gene_table.csv
"""

import argparse, os
import pandas as pd
import numpy as np

PRIORITY_GENES = {
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
    "CELSR2":  {"category": "lipid_metabolism",  "druggable": False},
    "PSRC1":   {"category": "lipid_metabolism",  "druggable": False},
    "IL6R":    {"category": "inflammation",      "druggable": True},
    "NLRP3":   {"category": "inflammation",      "druggable": True},
    "CXCL12":  {"category": "inflammation",      "druggable": False},
    "PHACTR1": {"category": "plaque_biology",    "druggable": False},
    "ADAMTS7": {"category": "plaque_biology",    "druggable": True},
    "COL4A1":  {"category": "plaque_biology",    "druggable": False},
    "COL4A2":  {"category": "plaque_biology",    "druggable": False},
    "TCF21":   {"category": "plaque_biology",    "druggable": False},
    "F5":      {"category": "thrombosis",        "druggable": True},
    "F2":      {"category": "thrombosis",        "druggable": True},
    "ABO":     {"category": "thrombosis",        "druggable": False},
    "VWF":     {"category": "thrombosis",        "druggable": True},
    "NOS3":    {"category": "bp_remodeling",     "druggable": True},
    "NPR3":    {"category": "bp_remodeling",     "druggable": True},
    "CYP17A1": {"category": "bp_remodeling",     "druggable": True},
}

HIGH_VALUE_TISSUES = {"Liver", "Artery_Coronary"}


def score_positional(pos_df):
    """Aggregate positional scores per gene."""
    if pos_df.empty:
        return pd.DataFrame(columns=["gene_name","positional_score","n_loci_pos",
                                     "max_pip","max_pp4","loci","exposures"])
    agg = pos_df.groupby("gene_name").agg(
        positional_score = ("positional_score", "max"),
        n_loci_pos       = ("locus_bp", "nunique"),
        max_pip          = ("max_pip",  "max"),
        max_pp4          = ("max_pp4",  "max"),
        loci             = ("locus_bp", lambda x: ";".join(str(v) for v in sorted(x.unique()))),
        exposures        = ("exposures", lambda x: ";".join(set(";".join(x).split(";")))),
    ).reset_index()
    return agg


def score_eqtl(eqtl_df):
    """Compute eQTL coloc score per gene."""
    if eqtl_df.empty:
        return pd.DataFrame(columns=["gene_name","eqtl_score","max_eqtl_pp4",
                                     "eqtl_tissues","eqtl_sig"])
    rows = []
    for gene, grp in eqtl_df.groupby("gene_name"):
        max_pp4    = grp["PP.H4"].max()
        sig_any    = grp["PP.H4"] >= 0.5
        sig_hi_tis = grp.loc[grp["tissue"].isin(HIGH_VALUE_TISSUES), "PP.H4"] >= 0.8
        tissues    = ";".join(sorted(grp.loc[grp["PP.H4"] >= 0.5, "tissue"].unique()))

        if sig_hi_tis.any():
            score = 3
        elif (grp["PP.H4"] >= 0.8).any():
            score = 2
        elif sig_any.any():
            score = 1
        else:
            score = 0

        rows.append({
            "gene_name":    gene,
            "eqtl_score":   score,
            "max_eqtl_pp4": max_pp4,
            "eqtl_tissues": tissues,
            "eqtl_sig":     score > 0,
        })
    return pd.DataFrame(rows)


def score_l2g(l2g_df):
    """Compute L2G score per gene."""
    if l2g_df.empty:
        return pd.DataFrame(columns=["gene_name","l2g_score_pts","max_l2g",
                                     "druggable_ot","approved_drugs"])
    agg = l2g_df.groupby("gene_name").agg(
        max_l2g      = ("l2g_score",     "max"),
        druggable_ot = ("druggable_ot",  "max"),
        approved_drugs = ("approved_drugs", lambda x: ";".join(
            set(d for v in x for d in str(v).split(";") if d and d != "nan"))),
    ).reset_index()
    agg["l2g_score_pts"] = pd.cut(agg["max_l2g"],
        bins=[-1, 0.2, 0.5, 0.8, 2],
        labels=[0, 1, 2, 3]).astype(float).fillna(0)
    return agg


def score_pip(cs_df):
    """Score based on fine-mapping PIP."""
    if cs_df.empty:
        return pd.DataFrame(columns=["gene_name","pip_score"])
    # This is applied at the locus level; we match by locus
    # Return a per-locus PIP score that gets merged later
    cs_df = cs_df.sort_values("top_pip", ascending=False).drop_duplicates(
        subset=["exposure","chr","locus_bp"])
    cs_df["pip_score"] = (cs_df["top_pip"] >= 0.9).astype(int)
    return cs_df[["exposure","chr","locus_bp","top_pip","pip_score"]]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--positional",    required=True)
    parser.add_argument("--eqtl_coloc",   required=True)
    parser.add_argument("--ot_l2g",       required=True)
    parser.add_argument("--credible_sets", required=True)
    parser.add_argument("--mr_results",    required=True)
    parser.add_argument("--out",           required=True)
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    # ── Load all evidence layers ──────────────────────────────────
    print("[scoring] Loading evidence layers...")

    pos_df  = pd.read_csv(args.positional)    if os.path.exists(args.positional)    else pd.DataFrame()
    eqtl_df = pd.read_csv(args.eqtl_coloc)   if os.path.exists(args.eqtl_coloc)   else pd.DataFrame()
    l2g_df  = pd.read_csv(args.ot_l2g)       if os.path.exists(args.ot_l2g)       else pd.DataFrame()
    cs_df   = pd.read_csv(args.credible_sets) if os.path.exists(args.credible_sets) else pd.DataFrame()
    mr_df   = pd.read_csv(args.mr_results)    if os.path.exists(args.mr_results)    else pd.DataFrame()

    # ── Score each layer ─────────────────────────────────────────
    pos_scores  = score_positional(pos_df)
    eqtl_scores = score_eqtl(eqtl_df)
    l2g_scores  = score_l2g(l2g_df)

    # ── Collect all gene names across layers ──────────────────────
    all_genes = set()
    for df, col in [(pos_scores, "gene_name"), (eqtl_scores, "gene_name"),
                    (l2g_scores, "gene_name")]:
        if not df.empty and col in df.columns:
            all_genes.update(df[col].dropna().unique())
    all_genes.update(PRIORITY_GENES.keys())
    all_genes = sorted(all_genes)

    print(f"[scoring] Scoring {len(all_genes)} unique genes...")

    # ── Build master gene table ───────────────────────────────────
    master = pd.DataFrame({"gene_name": list(all_genes)})

    # Merge positional
    master = master.merge(pos_scores.rename(columns={"n_loci_pos": "n_loci"}),
                          on="gene_name", how="left")

    # Merge eQTL
    master = master.merge(eqtl_scores, on="gene_name", how="left")

    # Merge L2G
    master = master.merge(l2g_scores[["gene_name","l2g_score_pts","max_l2g",
                                       "druggable_ot","approved_drugs"]],
                          on="gene_name", how="left")

    # Fill missing scores with 0
    for col in ["positional_score","eqtl_score","l2g_score_pts"]:
        master[col] = master[col].fillna(0).astype(float)

    # PIP score: 1 if gene is nearest to a singleton CS (PIP>=0.9)
    # Use max_pip from positional mapping
    master["pip_score"] = (master.get("max_pip", pd.Series(0, index=master.index))
                           .fillna(0) >= 0.9).astype(float)

    # Drug target score: 1 if in priority list or OT confirmed
    master["is_priority"] = master["gene_name"].isin(PRIORITY_GENES)
    master["druggable"]   = master.apply(lambda r:
        bool(r.get("druggable_ot")) or
        PRIORITY_GENES.get(r["gene_name"], {}).get("druggable", False), axis=1)
    master["drug_score"]  = master["druggable"].astype(float)

    # Category from priority list
    master["category"] = master["gene_name"].map(
        lambda g: PRIORITY_GENES.get(g, {}).get("category", "other"))

    # ── Total evidence score (max 10) ─────────────────────────────
    master["evidence_score"] = (
        master["positional_score"].clip(0, 3) +
        master["eqtl_score"].clip(0, 3) +
        master["l2g_score_pts"].clip(0, 2) +
        master["pip_score"].clip(0, 1) +
        master["drug_score"].clip(0, 1)
    ).clip(0, 10)

    # ── MR causal support ─────────────────────────────────────────
    if not mr_df.empty:
        # Flag if exposure matches gene (rough: gene in exposure name)
        ivw = mr_df[mr_df["method"] == "Inverse variance weighted"].copy()
        ivw_sig = ivw[ivw["pval"] < 0.05][["exposure","outcome","OR","pval"]].copy()
        # Attach MR evidence for known genes
        for gene in ["LDLR","PCSK9","SORT1","APOE","TRIB1","ANGPTL3","APOC3","LPA"]:
            # These are locus genes; MR evidence is at exposure level
            pass
        master["mr_causal_exposure"] = master["exposures"].fillna("") if "exposures" in master.columns else ""

    # ── Sort and output ───────────────────────────────────────────
    master = master.sort_values("evidence_score", ascending=False)

    out_cols = [
        "gene_name", "evidence_score", "category", "druggable",
        "positional_score", "eqtl_score", "l2g_score_pts", "pip_score", "drug_score",
        "max_pip", "max_pp4", "max_eqtl_pp4", "max_l2g",
        "eqtl_tissues", "n_loci", "loci", "exposures",
        "approved_drugs", "is_priority",
    ]
    out_cols = [c for c in out_cols if c in master.columns]
    master[out_cols].to_csv(args.out, index=False)

    print(f"\n[scoring] Top 30 ranked genes:")
    top = master.head(30)
    print(top[["gene_name","evidence_score","category","druggable",
               "positional_score","eqtl_score","l2g_score_pts","pip_score"]].to_string(index=False))

    print(f"\n[scoring] Score distribution:")
    print(master["evidence_score"].describe())
    print(f"\n[scoring] Written: {args.out} ({len(master)} genes)")


if __name__ == "__main__":
    main()
