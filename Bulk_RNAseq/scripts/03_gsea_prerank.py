#!/usr/bin/env python3
"""
Gene Set Encrichement Analysis of Bulk RNA-seq data
- Preranked GSEA (dual-species) using gseapy, with Enrichr GMT downloads.

Inputs:
  processed/Human_allgenes_stats.csv
  processed/Zebrafish_allgenes_stats.csv  
Outputs:
  processed/GSEA_human/
  processed/GSEA_zebrafish/

Notes:
  - Ranking metric: limma t-statistic (column 't')
  - Enrichr GMTs downloaded with UTF-8 to avoid encoding issues
  - Significance: FDR q-val < 0.25 (GSEA convention)
"""

from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import requests
import gseapy as gp


# ----------------------------
# Helpers
# ----------------------------
def build_rank_file(stats_path: Path, out_path: Path,
                    gene_col="gene", score_col="t") -> Path:
    """
    Build a prerank TSV: 2 columns (gene_name, score), no header.
    Dedup genes by max |score|, then sort descending.
    """
    df = pd.read_csv(stats_path)

    required = {gene_col, score_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {stats_path.name}: {missing}. "
                         f"Found: {list(df.columns)}")

    df = df.dropna(subset=[gene_col, score_col]).copy()
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=[score_col])

    # keep strongest |t| if duplicated gene symbols exist
    df = df.sort_values(score_col, key=lambda s: s.abs(), ascending=False)
    df = df.drop_duplicates(subset=gene_col, keep="first")

    rnk = df[[gene_col, score_col]].rename(columns={gene_col: "gene", score_col: "score"})
    rnk = rnk.sort_values("score", ascending=False)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    rnk.to_csv(out_path, sep="\t", index=False, header=False)

    print(f"[rank] wrote: {out_path} (n={len(rnk)})")
    return out_path


def ensure_rank_file(stats_path: Path, out_path: Path) -> Path:
    if out_path.exists():
        print(f"[rank] exists: {out_path} (skipping build)")
        return out_path
    return build_rank_file(stats_path, out_path)


def download_enrichr_gmt(library_name: str, out_path: Path) -> Path:
    """
    Download Enrichr gene set library as GMT (UTF-8).
    """
    url = "https://maayanlab.cloud/Enrichr/geneSetLibrary"
    params = {"mode": "text", "libraryName": library_name}
    r = requests.get(url, params=params, timeout=180)
    r.raise_for_status()
    out_path.write_text(r.text, encoding="utf-8")
    return out_path


def run_prerank(rnk_file: Path, out_dir: Path, libraries: list[str],
                permutation_num=2000, min_size=15, max_size=500,
                threads=4, seed=42) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    all_res = {}

    for gs in libraries:
        print("\n" + "=" * 80)
        print(f"[gmt] {gs}")
        print("=" * 80)

        gmt_path = out_dir / f"{gs}.gmt"
        download_enrichr_gmt(gs, gmt_path)
        print(f"[gmt] saved: {gmt_path}")

        print("\n" + "=" * 80)
        print(f"[gsea] prerank: {gs}")
        print("=" * 80)

        pre_res = gp.prerank(
            rnk=str(rnk_file),
            gene_sets=str(gmt_path),
            permutation_num=permutation_num,
            min_size=min_size,
            max_size=max_size,
            threads=threads,
            seed=seed,
            outdir=None,
            verbose=True,
        )

        res = pre_res.res2d.copy()
        for col in ["NES", "NOM p-val", "FDR q-val"]:
            if col in res.columns:
                res[col] = pd.to_numeric(res[col], errors="coerce")

        out_full = out_dir / f"GSEA_{gs}_full.csv"
        res.to_csv(out_full, index=False)

        sig = res[res["FDR q-val"] < 0.25].sort_values("FDR q-val")
        out_sig = out_dir / f"GSEA_{gs}_sig_FDR025.csv"
        sig.to_csv(out_sig, index=False)

        print(f"[save] full: {out_full}")
        print(f"[save] sig (FDR<0.25): {len(sig)} -> {out_sig}")

        all_res[gs] = res

    summary = pd.DataFrame([
        {
            "GeneSet": gs,
            "TotalTerms": len(df),
            "FDR<0.25": int((df["FDR q-val"] < 0.25).sum()),
            "FDR<0.10": int((df["FDR q-val"] < 0.10).sum()),
            "FDR<0.05": int((df["FDR q-val"] < 0.05).sum()),
        }
        for gs, df in all_res.items()
    ])
    return summary


# ----------------------------
# Analysis
# ----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--processed", type=Path, required=True,
                    help="Path to processed/ folder (contains *_allgenes_stats.csv)")
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()

    processed = args.processed

    # Human
    human_stats = processed / "Human_allgenes_stats.csv"
    human_out = processed / "GSEA_human"
    human_rnk = human_out / "Human_prerank_t.tsv"

    # Zebrafish
    z_stats = processed / "Zebrafish_allgenes_stats.csv"
    z_out = processed / "GSEA_zebrafish"
    z_rnk = z_out / "Zebrafish_prerank_t.tsv"

    # Libraries: 
    HUMAN_LIBS = [
        "MSigDB_Hallmark_2020",
        # "KEGG_2021_Human",
        # "Reactome_Pathways_2024",
        # "GO_Biological_Process_2025",
    ]
    ZFISH_LIBS = [
        "GO_Biological_Process_2018",
        "KEGG_2019",
        "WikiPathways_2018",
    ]

    print(f"[env] gseapy version: {gp.__version__}")

    if human_stats.exists():
        ensure_rank_file(human_stats, human_rnk)
        summ_h = run_prerank(human_rnk, human_out, HUMAN_LIBS, threads=args.threads)
        summ_h.to_csv(human_out / "GSEA_summary.csv", index=False)
        print("\n[human] summary\n", summ_h)
    else:
        print(f"[skip] missing: {human_stats}")

    if z_stats.exists():
        ensure_rank_file(z_stats, z_rnk)
        summ_z = run_prerank(z_rnk, z_out, ZFISH_LIBS, threads=args.threads)
        summ_z.to_csv(z_out / "GSEA_summary.csv", index=False)
        print("\n[zfish] summary\n", summ_z)
    else:
        print(f"[skip] missing: {z_stats}")


if __name__ == "__main__":
    main()