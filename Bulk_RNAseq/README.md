# Bulk RNA-seq analysis of SW480 and SW620 xenografts

This folder contains all scripts used to process and analyze bulk RNA-sequencing data from tumors dissected from 2dpi zebrafish xenografts of human colorectal cancer cell lines SW480 and SW620.

Dual-species transcript quantification was performed to separately analyze human tumor-derived transcripts and zebrafish host/tumor microenvironment (TME) transcripts.

---

## Overview of the analysis

1. **Transcript quantification**
   - Salmon was used to quantify transcripts independently against:
     - Human transcriptome (GRCh38)
     - Zebrafish transcriptome (GRCz11)
   - Lane-split FASTQ files were quantified as single-end reads.

2. **Gene-level summarization and differential expression**
   - Transcript-level abundances were summarized to gene-level counts using `tximport` with length-scaled TPMs.
   - Differential expression between SW620 and SW480 xenografts was performed using `edgeR` and `limma-voom`, with statistical testing via `limma-treat`.

3. **Gene set enrichment analysis (GSEA)**
   - Preranked GSEA was performed using gene-level t-statistics.
   - Enrichment was assessed separately for human tumor genes and zebrafish host genes.

---

## Scripts

Scripts should be run in the order listed below.

### 01_quantify_salmon_dual_species.sh
Quantifies RNA-seq reads using Salmon against human and zebrafish transcriptomes.

**Inputs**
- Lane-split FASTQ files (`.fastq.gz`)
- Reference transcriptomes (Ensembl cDNA FASTA files)

**Outputs**
- `results/salmon/human/<SAMPLE>/quant.sf`
- `results/salmon/zebrafish/<SAMPLE>/quant.sf`

---

### 02_DEG_analysis.R
Performs gene-level differential expression analysis for human and zebrafish data.

**Key steps**
- Transcript-to-gene mapping using Ensembl GTF files
- Gene-level count aggregation (`tximport`)
- Filtering of lowly expressed genes
- Normalization (TMM)
- Differential expression using `limma-voom` and `limma-treat`

**Outputs (written to `processed/`)**
- `<Species>_allgenes_stats.csv` (used for GSEA ranking)
- `<Species>_FULL_SW620_vs_SW480.csv`
- `<Species>_DEG_FDR0.05_LFC1_SW620_vs_SW480.csv`
- `<Species>_normalized_expression.csv`
- `processed_data.RData`
- `sessionInfo.txt`

---

### 03_GSEA_analysis.py
Performs preranked GSEA on human and zebrafish gene expression data.

**Key steps**
- Ranking genes by limma t-statistics
- Downloading gene sets from Enrichr
- Running preranked GSEA using `gseapy`
- Exporting full and significant enrichment results

---

## Software versions

- Salmon: v1.10.3
- R: v4.5.1
- edgeR: v4.8.2
- limma: v3.66.0
- tximport: v1.38.2
- gseapy: see `sessionInfo.txt` / Python environment

---

## Data availability

Raw RNA-seq data is deposited in GEO under accession **GSE163746**.

Processed quantification files (`quant.sf`) and downstream analysis outputs generated for this study are provided as Supplementary Data accompanying the manuscript.

This repository contains only the scripts required to reproduce the analyses.