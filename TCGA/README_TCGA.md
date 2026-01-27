
# TCGA survival analysis of CD24 in colorectal cancer

This folder contains scripts used to analyze publicly available TCGA colorectal cancer datasets (COAD and READ) to evaluate the prognostic value of CD24 expression, alone and in combination with a SIGLEC10-associated tumor-associated macrophage (TAM) signature.

---

## Overview of the analysis

### 1. Data acquisition

Gene expression, clinical, and survival data were downloaded from the UCSC Xena Browser for:
- Colon Adenocarcinoma (COAD)
- Rectum Adenocarcinoma (READ)

### 2. Patient-level cohort construction

Expression and clinical datasets were harmonized and merged prior to analysis.
- For patients with multiple samples a single representative sample was selected to ensure statistical independence.
- Only patients with valid survival and staging information were retained.
- Analyses were restricted to AJCC stages I–IV.

### 3. Expression stratification

#### CD24 alone
- CD24 expression was extracted from normalized RNA-seq expression matrices.
- Patients were stratified into **CD24 High** and **CD24 Low** groups using a cohort-wide median split.

#### CD24 + SIGLEC10+ TAM signature
- A macrophage score was computed using pan-macrophage marker genes.
- A SIGLEC10+ TAM score was calculated as: SIGLEC10 expression × macrophage score
- Patients were stratified into four groups using median splits:
    - Low CD24 / Low SIGLEC10+ TAM
    - Low CD24 / High SIGLEC10+ TAM
    - High CD24 / Low SIGLEC10+ TAM
    - High CD24 / High SIGLEC10+ TAM

### 4. Survival analysis

- Overall Survival (OS) and Progression-Free Interval (PFI) were analyzed using Kaplan–Meier estimates.
- All survival analyses were **truncated at 5 years** to reduce bias from sparse long-term follow-up.
- Statistical significance was assessed using:
    - Log-rank tests (CD24 High vs Low)
    - Multivariate log-rank tests (CD24/SIGLEC10 4-group analysis)
- Cox proportional hazards models were fitted where possible, but unstable or non-estimable hazard ratios (e.g. due to low event counts) were omitted.

---
## Scripts

Scripts should be run in the order listed below:

### (1) `tcga_analysis.py`

**Inputs** - Downloaded from UCSC Xena Browswer (https://xenabrowser.net/).

Gene expression matrices 
- COAD: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-COAD.star_fpkm-uq.tsv.gz
- READ: https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-READ.star_fpkm-uq.tsv.gz

PanCancer survival data: https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/Survival_SupplementalTable_S1_20171025_xena_sp
(from all TCGA datasets - filtered for COAD+READ during analysis)

**Purpose**
- Load and harmonize TCGA COAD and READ datasets
- Perform patient-level sample merging
- Compute CD24 and SIGLEC10+ TAM scores
- Output a single analysis-ready table

**Key steps**
- Expression and clinical data loading
- Patient-level aggregation
- Stage filtering (I–IV)
- SIGLEC10+ TAM score computation
- Output of a unified analysis table

**Output**
- `results/processed/tcga_analysis_table.csv`


### (2) `tcga_figures.py`

**Purpose**
- Generate Kaplan–Meier survival figures
- Create a master survival summary table
- Produce a comprehensive clinical characteristics table

**Key outputs**
Written to `results/figures_tables/`:

- `CD24_OS_PFI_5yr.png`  
Stage-stratified OS and PFI Kaplan–Meier plots for CD24 High vs Low

- `CD24_SIGLEC10_OS_PFI_5yr.png`  
Stage-stratified OS and PFI Kaplan–Meier plots for CD24/SIGLEC10+ TAM groups

- `Master_Survival_Summary.csv`  
Single consolidated survival summary table including:
- Overall cohort statistics
- Stage-stratified log-rank p-values
- Group sizes and event counts
- Hazard ratios (when estimable)

- `AllPatients_ClinicalCharacteristics.png`  
- `AllPatients_ClinicalCharacteristics.csv`  

---
## Data availability

All TCGA data used in this analysis are publicly available through the UCSC Xena Browser.
Processed outputs are provided as Supplementary Data accompanying the manuscript.
This repository contains only the scripts required to reproduce the analyses.
Software and package versions used for the analysis are documented in the project environment file.