# Zebrahub macrophage analysis (zebrafish atlas)

This folder contains a lightweight analysis of publicly available single-cell RNA-seq
zebrafish atlases (ZebraHub) to examine expression of **cd24a-like / siglec-related genes**
in macrophages across early developmental stages.

The goal of this analysis is **illustrative and supportive**, not to generate a new atlas:
we use existing curated datasets to show that macrophages express relevant immune
regulatory genes during early zebrafish development, supporting observations made
in tumor xenograft experiments.

---

## Input data (ZebraHub atlases)

The analysis uses preprocessed `.h5ad` files provided by **ZebraHub**:

- **2 dpf atlas**
- **3 dpf atlas**
- **5 dpf atlas**

You can download them from:

🔗 https://zebrahub.ds.czbiohub.org  

Navigate to **Downloads → Cell atlases → scRNA-seq atlases**, and download the
corresponding `.h5ad` files for:
- `2 dpf`
- `3 dpf`
- `5 dpf`

These files already contain:
- Normalized counts
- Cell type / cluster annotations
- Developmental stage metadata

No raw FASTQs are required.

---

## Script

### `zebrahub_analysis_figures.py`

This script reproduces all analyses and figures previously generated in the notebook.

**Main steps:**
1. Load 2, 3, and 5 dpf `.h5ad` atlases
2. Concatenate datasets and harmonize metadata
3. Identify macrophages using canonical markers:
   - `mpeg1.1`
   - `csf1ra`
   - `csf1rb`
   - `mfap4`
4. Normalize and log-transform expression
5. Visualize candidate genes using:
   - Violin plots
   - Dot plots
6. Export a summary table with mean expression per stage

---

## How to run

```bash
python zebrahub_analysis_figures.py \
  --h5ad-2dpf path/to/2dpf_atlas.h5ad \
  --h5ad-3dpf path/to/3dpf_atlas.h5ad \
  --h5ad-5dpf path/to/5dpf_atlas.h5ad \
  --outdir results/zebrahub
```

### Optional arguments
- `--mac-markers`  
  Comma-separated list of macrophage marker genes  
  *(default: mpeg1.1,csf1ra,csf1rb,mfap4)*

- `--mac-threshold`  
  Mean expression threshold used to label macrophages  
  *(default: 2.0)*

- `--candidate-genes`  
  Genes to visualize (comma-separated)

---

## Outputs

Written to the output directory:

- **Violin plots** of candidate gene expression in macrophages
- **Dot plots** summarizing expression across stages
- **CSV table** with mean expression per gene and developmental stage

Example outputs:
```
results/zebrahub/
├── Violin_candidate_genes.png
├── Dotplot_candidate_genes.png
└── Mean_expression_macrophages.csv
```

---

## Notes on interpretation

- This analysis is **descriptive** and leverages curated cell atlases.
- Macrophages are identified computationally using marker expression,
  not re-clustered.
- Results support the presence of immune-regulatory gene expression
  in zebrafish macrophages during early development.
- These data are intended to complement xenograft and bulk RNA-seq analyses,
  not to serve as an independent discovery dataset.

---

## Citation

If you use the ZebraHub atlases, please cite the original ZebraHub publication
and data resource as indicated on their website.

🔗 https://zebrahub.ds.czbiohub.org
