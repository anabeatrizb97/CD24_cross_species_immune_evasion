# TCGA survival analysis of CD24 in colorectal cancer

This folder contains scripts used to analyze publicly available TCGA colorectal cancer datasets (COAD and READ) to assess the prognostic value of CD24 expression, alone and in combination with SIGLEC10-associated tumor-associated macrophage (TAM) signatures.
All analyses were performed using processed, publicly available TCGA data and did not involve raw sequencing files.

---

## Overview of the analysis

1. **Data acquisition**

Gene expression, clinical and survival data were downloaded from the UCSC Xena Browser for:
- Colon Adenocarcinoma (COAD)
- Rectum Adenocarcinoma (READ)

2. **Expression stratification**

- CD24 only:
CD24 expression was extracted from normalized RNA-seq expression matrices. Samples were stratified into high and low expression groups based on cohort-specific thresholds (median split).

- CD24 + Siglec-10 combined:
Combined stratification was performed to assess the interaction between CD24 expression and macrophage-expressed Siglec-10. SIGLEC10+ TAM score was computed using a predefined macrophage-related gene set.

3. **Survival analysis**

Overall survival was analyzed using Kaplan–Meier estimates and statistical significance was assessed using log-rank tests.

---

## Script

This script should only be ran after downloading the input files bellow.
### tcga_analysis.py

# Key steps

- Loading and harmonizing expression and clinical data
- Sample filtering and cohort merging (COAD + READ)
- SIGLEC10+ TAM score computation
- Stratification by markers' expression
- Kaplan–Meier survival analysis and log-rank testing

---

## Inputs 

Downloaded from UCSC Xena Browswer (https://xenabrowser.net/).

Gene expression matrices 
- COAD: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.star_fpkm-uq.tsv.gz
- READ: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-READ.star_fpkm-uq.tsv.gz
(at the time of analysis, COAD n=514 and READ n=177)

Clinical data:
- COAD: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.clinical.tsv.gz
- READ: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-READ.clinical.tsv.gz
(at the time of analysis, COAD n=562 and READ n=187)

PanCancer survival data: https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp
(from all TCGA datasets - filtered for COAD+READ during analysis)

---

## Outputs

- Survival analyses summary tables
- Kaplan–Meier plots
- Detailed patient data table

---

## Software versions

- Python:3.11
- pandas: v2.1.0
- numpy: v2.3.4
- lifelines: 0.30.0
- matplotlib: v3.8.0

---

## Data availability

All TCGA data used in this analysis are publicly available through the UCSC Xena Browser.
Processed outputs and figures are provided as Supplementary Data accompanying the manuscript.
This repository contains only the scripts required to reproduce the analyses.