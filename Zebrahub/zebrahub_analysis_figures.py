#!/usr/bin/env python3
"""
ZebraHub zebrafish atlas scRNAseq analysis
- Gene expression analysis at specific developmental timepoints

Analysis overview
------------
1) Loads Zebrahub 2dpf, 3dpf, 5dpf h5ad gene expression atlases.
2) Detects macrophage-enriched (cluster, timepoint) pairs based on expression of macrophage markers.
3) Creates figures:
   - Violin plot of the expression of candidate CD24 receptor genes in zebrafish macrophages
   - DotPlot across selected cluster × timepoint groups for candidate genes + macrophage markers (and csv table)

Inputs
------
Zebrahub 2dpf, 3dpf, 5dpf .h5ad files under:

<project-root>/data/
- `zf_atlas_*2dpf*_release*.h5ad`
- `zf_atlas_*3dpf*_release*.h5ad`
- `zf_atlas_*5dpf*_release*.h5ad`

Outputs
-------
Written to:
  <project-root>/results/

"""

from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------
# Helpers
# -----------------------------

def load_and_concat(h5ads: list[Path]) -> sc.AnnData:
    adatas = []
    for fp in h5ads:
        print(f"[load] {fp}")
        ad = sc.read_h5ad(fp)

        if "timepoint_cluster" not in ad.obs.columns:
            raise KeyError(f"{fp.name}: missing obs['timepoint_cluster']")
        if "timepoint" not in ad.obs.columns:
            raise KeyError(f"{fp.name}: missing obs['timepoint']")

        ad.obs["timepoint_cluster"] = ad.obs["timepoint_cluster"].astype(str)
        ad.obs["timepoint"] = ad.obs["timepoint"].astype(str)

        adatas.append(ad)

    ad_combined = sc.concat(adatas, join="inner", merge="same", label="source", keys=[p.stem for p in h5ads])
    print(f"[concat] combined shape: {ad_combined.shape}")
    return ad_combined


def subset_macrophages(ad_combined: sc.AnnData) -> sc.AnnData:
    print("Subsetting macrophages...")

    mac_markers = ["mpeg1.1", "csf1ra", "csf1rb", "mfap4"]
    present = [g for g in mac_markers if g in ad_combined.var_names]
    missing = [g for g in mac_markers if g not in ad_combined.var_names]

    if missing:
        print(f"[warn] missing macrophage markers in var_names (ignoring): {missing}")
    if not present:
        raise ValueError("None of the macrophage markers were found in var_names; cannot subset macrophages.")

    threshold = 2.0 

    df_markers = sc.get.obs_df(ad_combined, keys=present)
    df_markers["timepoint_cluster"] = ad_combined.obs["timepoint_cluster"].astype(str).values
    df_markers["timepoint"] = ad_combined.obs["timepoint"].astype(str).values

    # Mean expression per (cluster, timepoint)
    cluster_expr = df_markers.groupby(["timepoint_cluster", "timepoint"], observed=True)[present].mean()

    mac_clusters_auto_tp = cluster_expr[(cluster_expr > threshold).any(axis=1)].index.tolist()
    print("\nAutomatically detected macrophage clusters (cluster, timepoint):")
    print(mac_clusters_auto_tp)

    if not mac_clusters_auto_tp:
        raise ValueError(
            "No macrophage clusters detected with the current threshold.\n"
            "Try lowering threshold (e.g. 1.0) or confirm markers exist / are normalized as expected."
        )

    # Subset cells matching selected (cluster,timepoint)
    obs_pairs = list(zip(ad_combined.obs["timepoint_cluster"].astype(str), ad_combined.obs["timepoint"].astype(str)))
    keep = np.array([p in set(mac_clusters_auto_tp) for p in obs_pairs], dtype=bool)

    macrophages = ad_combined[keep].copy()
    print(f"\nMacrophage subset shape: {macrophages.shape}")
    print("Timepoints included:", sorted(macrophages.obs["timepoint"].unique()))

    return macrophages


def normalize_log1p(adata: sc.AnnData) -> None:
    print("Normalizing and log-transforming ...")
    if "log1p" not in adata.uns:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    print("Normalization and log-transform complete.")


def make_violin(macrophages: sc.AnnData, out_png: Path) -> None:
    genes = ["siglec15l", "mag", "si:dkey-24p1.7"]
    genes_in_data = [g for g in genes if g in macrophages.var_names]
    print(f"Genes found in dataset (violin): {genes_in_data}")
    if not genes_in_data:
        raise ValueError("None of the candidate genes were found in var_names; cannot make violin plot.")

    fig, ax = plt.subplots(figsize=(4, 6))  
    sc.pl.violin(
        macrophages,
        keys=genes,  
        rotation=45,
        stripplot=True,
        jitter=0.4,
        ax=ax,
        show=False,
        palette=["#F4C2C2", "#B0C4D9", "#D67373"],  
    )

    ax.set_ylabel("")
    ax.text(
        -0.20,
        0.58,
        "Gene expression in macrophages",
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="center",
        ha="right",
        rotation=90,
    )
    ax.text(
        -0.15,
        0.58,
        "(log-normalized)",
        transform=ax.transAxes,
        fontsize=10,
        fontweight="normal",
        va="center",
        ha="right",
        rotation=90,
    )

    ax.set_xlabel("")
    ax.tick_params(axis="both", labelsize=9)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["left"].set_visible(True)
    ax.grid(False)

    plt.tight_layout()
    fig.canvas.draw()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[save] {out_png}")


def make_dotplot(ad_combined: sc.AnnData, out_png: Path) -> None:
    genes_of_interest = ["siglec15l", "mag", "si:dkey-24p1.7", "mpeg1.1", "csf1ra", "mfap4"]
    var_genes = ["siglec15l", "mag", "si:dkey-24p1.7"]

    # Ensure required obs fields exist
    if "timepoint_cluster" not in ad_combined.obs.columns or "timepoint" not in ad_combined.obs.columns:
        raise KeyError("ad_combined.obs must contain 'timepoint_cluster' and 'timepoint'.")

    # Build cluster_timepoint in-place
    ad_combined.obs["timepoint_cluster"] = ad_combined.obs["timepoint_cluster"].astype(str)
    ad_combined.obs["timepoint"] = ad_combined.obs["timepoint"].astype(str)
    ad_combined.obs["cluster_timepoint"] = ad_combined.obs["timepoint_cluster"] + "_" + ad_combined.obs["timepoint"]

    # Only keep genes present
    present_var = [g for g in var_genes if g in ad_combined.var_names]
    if not present_var:
        raise ValueError("None of the dotplot variable genes were found in var_names; cannot make dotplot.")

    present_interest = [g for g in genes_of_interest if g in ad_combined.var_names]
    if not present_interest:
        raise ValueError("None of the dotplot genes were found in var_names; cannot make dotplot.")

    df_small = sc.get.obs_df(ad_combined, keys=list(set(present_interest + present_var)))
    df_small["cluster_timepoint"] = ad_combined.obs["cluster_timepoint"].values

    # Mean expression per cluster_timepoint for var genes
    cluster_means = df_small.groupby("cluster_timepoint", observed=True)[present_var].mean()
    cluster_vars = cluster_means.var(axis=1)

    if "siglec15l" in df_small.columns:
        siglec_means = df_small.groupby("cluster_timepoint", observed=True)["siglec15l"].mean()
        siglec_clusters = siglec_means[siglec_means > 0.001].index.tolist()
    else:
        siglec_clusters = []

    # Top 10 variable clusters
    top_clusters = cluster_vars.sort_values(ascending=False).head(10).index.tolist()

    for c in siglec_clusters[:2]:
        if c not in top_clusters:
            top_clusters.append(c)

    # -----------------------------
    # Macrophage clusters
    # -----------------------------
    mac_clusters_tp = [("22", "2dpf"), ("39", "3dpf"), ("17", "5dpf")]
    mac_clusters_tp_str = [f"{c}_{t}" for c, t in mac_clusters_tp]

    # Sort clusters
    mac_first_sorted = [
        c
        for tp in ["2dpf", "3dpf", "5dpf"]
        for c in mac_clusters_tp_str
        if c.endswith(tp) and c in top_clusters
    ]
    siglec_sorted = [c for c in siglec_clusters if c in top_clusters and c not in mac_first_sorted]
    other_clusters_sorted = [c for c in top_clusters if c not in mac_first_sorted + siglec_sorted]
    all_clusters = mac_first_sorted + siglec_sorted + other_clusters_sorted

    mask = ad_combined.obs["cluster_timepoint"].isin(top_clusters)
    adata_plot = ad_combined[mask].copy()

    existing = list(pd.unique(adata_plot.obs["cluster_timepoint"]))
    categories_order = [c for c in all_clusters if c in existing]

    dp = sc.pl.dotplot(
        adata_plot,
        var_names=present_interest,
        groupby="cluster_timepoint",
        standard_scale="var",
        swap_axes=True,
        dot_max=0.9,
        categories_order=categories_order,
        show=False,
        return_fig=True,
    )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    dp.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"[save] {out_png}")


# -----------------------------
# Main
# -----------------------------

def main():
    script_dir = Path(__file__).resolve().parent

    data_dir = script_dir / "data"
    out_dir = script_dir / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Detect inputs
    h5ads = sorted(data_dir.glob("zf_atlas_*dpf*_release.h5ad"))
    if not h5ads:
        raise FileNotFoundError(
            f"No .h5ad inputs found.\n"
            f"Looked in: {data_dir}\n"
            f"Expected pattern: zf_atlas_*dpf_v1_release.h5ad\n"
            f"Example: zf_atlas_2dpf_v1_release.h5ad"
        )

    print(f"[paths] script_dir: {script_dir}")
    print(f"[paths] data_dir:   {data_dir}")
    print(f"[paths] out_dir:    {out_dir}")
    print("[inputs]")
    for f in h5ads:
        print(" -", f.name)

    ad_combined = load_and_concat(h5ads)

    macrophages = subset_macrophages(ad_combined)

    print("Defining candidate receptor genes ...")
    candidate_genes = ["mag", "siglec15l", "si:dkey-24p1.7"]
    genes_in_data = [g for g in candidate_genes if g in macrophages.var_names]
    print(f"Genes found in dataset: {genes_in_data}")

    normalize_log1p(macrophages)

    make_violin(macrophages, out_dir / "candidates_macrophageexpression_violinplot.png")
    make_dotplot(ad_combined, out_dir / "candidate_genes_dotplot.png")

    print("\nDone.")


if __name__ == "__main__":
    main()