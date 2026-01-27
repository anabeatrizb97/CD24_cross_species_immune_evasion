#!/usr/bin/env python3
"""
Gene Set Encrichement Analysis (GSEA) of dual-species bulk RNA-seq data
- Preranked GSEA using gseapy, with Enrichr API for genesets.

Inputs:
  processed/Human_allgenes_stats.csv
  processed/Zebrafish_allgenes_stats.csv  
Outputs:
  processed/GSEA_human/
  processed/GSEA_zebrafish/
  refs/gsea_snapshots/  (GMT snapshots of the *API-fetched* libraries)

Key settings:
  - Ranking metric: limma t-statistic (column 't')
  - Permutations: 2000
  - Gene set size: 15–500
  - Significance: FDR q-val < 0.25

"""

from __future__ import annotations

from pathlib import Path
import re
import sys

import pandas as pd
import gseapy as gp


# ----------------------------
# Paths
# ----------------------------
PROJECT_ROOT = Path(__file__).resolve().parent 
PROCESSED_DIR = PROJECT_ROOT / "processed"

# Snapshots (API -> GMT written here)
SNAPSHOT_DIR = PROJECT_ROOT / "refs" / "gsea_snapshots"
SNAPSHOT_HUMAN = SNAPSHOT_DIR / "human"
SNAPSHOT_FISH = SNAPSHOT_DIR / "zebrafish"

# Inputs
HUMAN_STATS = PROCESSED_DIR / "Human_allgenes_stats.csv"
FISH_STATS = PROCESSED_DIR / "Zebrafish_allgenes_stats.csv"

# Outputs
OUT_HUMAN = PROCESSED_DIR / "GSEA_human"
OUT_FISH = PROCESSED_DIR / "GSEA_zebrafish"

# Rank files
RNK_HUMAN = OUT_HUMAN / "Human_prerank_t.tsv"
RNK_FISH = OUT_FISH / "Zebrafish_prerank_t.tsv"

# Parameters (match what you’ve been using)
PERMUTATIONS = 2000
MIN_SIZE = 15
MAX_SIZE = 500
THREADS = 4
SEED = 42

FDR_SIG = 0.25


# ----------------------------
# Helpers
# ----------------------------
def build_rank_file(stats_csv: Path, out_tsv: Path, gene_col: str = "gene", score_col: str = "t") -> None:
    """
    Build prerank TSV: 2 columns, no header: gene, score.
    Deduplicate by keeping max |score| per gene; then sort descending.
    """
    df = pd.read_csv(stats_csv)
    missing = {gene_col, score_col} - set(df.columns)
    if missing:
        raise ValueError(f"{stats_csv.name} missing columns {missing}. Found: {list(df.columns)}")

    df = df.dropna(subset=[gene_col, score_col]).copy()
    df[score_col] = pd.to_numeric(df[score_col], errors="coerce")
    df = df.dropna(subset=[score_col])

    df = df.sort_values(score_col, key=lambda s: s.abs(), ascending=False)
    df = df.drop_duplicates(subset=gene_col, keep="first")

    rnk = df[[gene_col, score_col]].rename(columns={gene_col: "gene", score_col: "score"})
    rnk = rnk.sort_values("score", ascending=False)

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    rnk.to_csv(out_tsv, sep="\t", index=False, header=False)
    print(f"[rank] wrote: {out_tsv} (n={len(rnk)})")


def read_gene_list(rnk_file: Path) -> list[str]:
    rnk_df = pd.read_csv(rnk_file, sep="\t", header=None, names=["gene", "score"])
    return rnk_df["gene"].astype(str).tolist()


def write_gmt_snapshot(gmt_dict: dict[str, list[str]], out_gmt: Path) -> None:
    out_gmt.parent.mkdir(parents=True, exist_ok=True)
    with out_gmt.open("w", encoding="utf-8") as f:
        for term, genes in gmt_dict.items():
            f.write(term + "\tNA\t" + "\t".join(map(str, genes)) + "\n")
    print(f"[snapshot] saved GMT: {out_gmt}")


def run_prerank_api(
    *,
    species_name: str,
    organism: str,
    rnk_file: Path,
    libraries: list[str],
    out_dir: Path,
    snapshot_dir: Path,
) -> None:

    print("\n" + "=" * 80)
    print(f"[{species_name}] starting (organism={organism})")
    print("=" * 80)

    gene_list = read_gene_list(rnk_file)

    for lib in libraries:
        print("\n" + "-" * 80)
        print(f"[{species_name}] library: {lib}")
        print("-" * 80)

        gmt_dict = gp.get_library(
            name=lib,
            organism=organism,
            min_size=MIN_SIZE,
            max_size=MAX_SIZE,
            gene_list=gene_list,
        )
        print(f"[api] fetched gene sets after filtering: {len(gmt_dict)}")

        # Geneset gmt snapshot
        snapshot_path = snapshot_dir / f"{lib}.{organism}.min{MIN_SIZE}_max{MAX_SIZE}.gmt"
        write_gmt_snapshot(gmt_dict, snapshot_path)

        # Run prerank
        pre_res = gp.prerank(
            rnk=str(rnk_file),
            gene_sets=gmt_dict,
            organism=organism,
            permutation_num=PERMUTATIONS,
            min_size=MIN_SIZE,
            max_size=MAX_SIZE,
            threads=THREADS,
            seed=SEED,
            outdir=None,
            verbose=True,
        )

        res = pre_res.res2d.copy()
        for col in ["NES", "NOM p-val", "FDR q-val"]:
            if col in res.columns:
                res[col] = pd.to_numeric(res[col], errors="coerce")

        out_dir.mkdir(parents=True, exist_ok=True)

        full_csv = out_dir / f"GSEA_{lib}_full.csv"
        sig_csv = out_dir / f"GSEA_{lib}_sig_FDR{FDR_SIG}.csv"

        res.to_csv(full_csv, index=False)
        sig = res[res["FDR q-val"] < FDR_SIG].sort_values("FDR q-val")
        sig.to_csv(sig_csv, index=False)

        print(f"[save] full: {full_csv}")
        print(f"[save] sig (FDR<{FDR_SIG}): {len(sig)} -> {sig_csv}")

    print(f"\n[{species_name}] done.")


# ----------------------------
# Main
# ----------------------------
def main() -> None:
    print(f"[env] gseapy version: {gp.__version__}")
    print(f"[paths] PROJECT_ROOT: {PROJECT_ROOT}")
    print(f"[paths] PROCESSED_DIR: {PROCESSED_DIR}")
    print(f"[paths] SNAPSHOT_DIR: {SNAPSHOT_DIR}")

    for p in [HUMAN_STATS, FISH_STATS]:
        if not p.exists():
            sys.exit(f"ERROR: missing input file: {p}")

    # Build rank files if needed
    if not RNK_HUMAN.exists():
        build_rank_file(HUMAN_STATS, RNK_HUMAN)
    else:
        print(f"[rank] exists: {RNK_HUMAN} (skipping)")

    if not RNK_FISH.exists():
        build_rank_file(FISH_STATS, RNK_FISH)
    else:
        print(f"[rank] exists: {RNK_FISH} (skipping)")

    HUMAN_LIBS = [
        "MSigDB_Hallmark_2020",
    ]

    FISH_LIBS = [
        "GO_Biological_Process_2018",
        "KEGG_2019",
        "WikiPathways_2018",
    ]

    # Run
    run_prerank_api(
        species_name="Human",
        organism="Human",
        rnk_file=RNK_HUMAN,
        libraries=HUMAN_LIBS,
        out_dir=OUT_HUMAN,
        snapshot_dir=SNAPSHOT_HUMAN,
    )

    run_prerank_api(
        species_name="Zebrafish",
        organism="Fish",
        rnk_file=RNK_FISH,
        libraries=FISH_LIBS,
        out_dir=OUT_FISH,
        snapshot_dir=SNAPSHOT_FISH,
    )

    print("\nAll GSEA finished.")


if __name__ == "__main__":
    main()