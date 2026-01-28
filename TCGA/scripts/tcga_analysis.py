#!/usr/bin/env python3
"""
TCGA COAD/READ survival analysis 
- Acording to CD24 expression 

This script builds a single, harmonized analysis table used for all TCGA-based survival analyses in the study. It integrates gene expression 
data with clinical and survival annotations and performs all preprocessing required prior to statistical testing and figure generation.

Inputs (from USC Xena Browser):
  (1) COAD and READ patient-level gene expression matrices (`x.star_fpkm-uq.tsv`; genes x samples)
  (2) PanCancer survival table (patient-level clinical + survival data)

Output:
  (1) Merged analysis table (for downstream analysis with `tcga_figures.py`)
  Contains patient-level:
  • Clinical annotations
  • Overall Survival (OS) and Progression-Free Interval (PFI)
  • CD24 expression

Path structure:
<PROJECT_ROOT>/
  tcga_analysis.py
  tcga_figure.py
  data/
    TCGA-COAD.star_fpkm-uq.tsv
    TCGA-READ.star_fpkm-uq.tsv
    TCGA_PanCancer_Survival.txt
  results/
    processed/tcga_analysis_table.csv

"""
from __future__ import annotations
from pathlib import Path
import re
import numpy as np
import pandas as pd

DAYS_PER_MONTH = 30.4375

# 5-year cutoff (in months) for PFI censoring columns
PFI_CUTOFF_MONTHS = 60.0

# CD24 Ensembl ID
GENES = {
    "ENSG00000272398": "CD24",
  }

_STAGE_ORDER = {"I": 1, "II": 2, "III": 3, "IV": 4}


def project_root() -> Path:
    return Path(__file__).resolve().parent


def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)


def tcga_patient_id(barcode: str) -> str:
    if pd.isna(barcode):
        return np.nan
    return str(barcode)[:12]


def tcga_sample_type(barcode: str) -> str:
    """
    TCGA sample type is the first 2 digits of the 4th barcode field.
    Example: TCGA-XX-YYYY-01A-... -> sample_type = "01" (Primary Tumor)
    """
    if pd.isna(barcode):
        return np.nan
    parts = str(barcode).split("-")
    if len(parts) >= 4 and len(parts[3]) >= 2:
        return parts[3][:2]
    return np.nan


def normalize_stage(raw):
    if pd.isna(raw):
        return np.nan
    s = str(raw).upper().strip()
    m = re.search(r"STAGE\s*([IV]+)", s)
    if m:
        return m.group(1)
    if s in {"I", "II", "III", "IV"}:
        return s
    return np.nan


def pick_max_stage(x: pd.Series):
    x = x.dropna().astype(str)
    x = x[x.isin(_STAGE_ORDER)]
    if len(x) == 0:
        return np.nan
    return max(x, key=lambda s: _STAGE_ORDER.get(s, 0))


def first_nonnull(x: pd.Series):
    x = x.dropna()
    return x.iloc[0] if len(x) else np.nan


def _base_ensembl_id(idx: pd.Index) -> pd.Index:
    s = pd.Index(idx.astype(str))
    s = s.str.split("|").str[0]
    s = s.str.split(".").str[0]
    return s


def load_xena_expression_minimal(path: Path, needed_ensg: list[str]) -> pd.DataFrame:
    """
    Xena expression matrix is genes x samples.
    We normalize gene IDs, subset to needed genes, then transpose to samples x genes.
    """
    raw = pd.read_csv(path, sep="\t", index_col=0)
    raw.index = _base_ensembl_id(raw.index)

    # If stripping versions creates duplicates, average them
    if raw.index.has_duplicates:
        raw = raw.groupby(raw.index).mean(numeric_only=True)

    raw = raw.reindex(needed_ensg)

    df = raw.T
    df.index.name = "sample"
    df.reset_index(inplace=True)
    return df


def load_survival(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df.rename(columns={df.columns[0]: "sample"}, inplace=True)

    df["patient_id"] = df["sample"].map(tcga_patient_id)

    df["time"] = pd.to_numeric(df.get("OS.time"), errors="coerce") / DAYS_PER_MONTH
    df["event"] = pd.to_numeric(df.get("OS"), errors="coerce")

    df["pfi_time"] = pd.to_numeric(df.get("PFI.time"), errors="coerce") / DAYS_PER_MONTH
    df["pfi_event"] = pd.to_numeric(df.get("PFI"), errors="coerce")

    if "ajcc_pathologic_tumor_stage" in df.columns:
        df["stage"] = df["ajcc_pathologic_tumor_stage"]
    elif "stage" in df.columns:
        df["stage"] = df["stage"]
    else:
        df["stage"] = np.nan

    df["stage"] = df["stage"].map(normalize_stage)

    wanted = [
        "sample",
        "patient_id",
        "stage",
        "time",
        "event",
        "pfi_time",
        "pfi_event",
        "gender",
        "race",
        "vital_status",
        "age_at_initial_pathologic_diagnosis",
        "cancer type abbreviation",
    ]
    for c in wanted:
        if c not in df.columns:
            df[c] = np.nan

    return df[wanted]


def pick_representative_sample(sub: pd.DataFrame, gene_cols: list[str]) -> pd.Series:
    """
    Choose one main sample per patient:
      1) Prefer Primary Tumor (sample_type == "01") if present
      2) Within the remaining, choose the sample with highest CD24 expression
         (falls back to highest median across loaded genes if CD24 missing)
    """
    sub = sub.copy()

    # Prefer primary tumor if any
    if (sub["sample_type"] == "01").any():
        sub = sub[sub["sample_type"] == "01"].copy()

    cd24_ensg = "ENSG00000272398"
    if cd24_ensg in sub.columns:
        key = pd.to_numeric(sub[cd24_ensg], errors="coerce")
    else:
        key = sub[gene_cols].median(axis=1, numeric_only=True)

    # If all NaN, fall back deterministically to first row
    if key.notna().any():
        return sub.loc[key.idxmax()]
    return sub.iloc[0]


def main() -> None:
    root = project_root()
    data_dir = root / "data"
    out_dir = root / "results" / "processed"
    ensure_dir(out_dir)

    needed_ensg = list(GENES.keys())

    print("Loading expression matrices...")
    expr_coad = load_xena_expression_minimal(data_dir / "TCGA-COAD.star_fpkm-uq.tsv", needed_ensg)
    expr_read = load_xena_expression_minimal(data_dir / "TCGA-READ.star_fpkm-uq.tsv", needed_ensg)
    expr = pd.concat([expr_coad, expr_read], ignore_index=True)

    expr["patient_id"] = expr["sample"].map(tcga_patient_id)
    expr["sample_type"] = expr["sample"].map(tcga_sample_type)

    gene_cols = [c for c in expr.columns if c not in {"sample", "patient_id", "sample_type"}]
    print(f"  Expression rows (samples): {expr.shape[0]} | genes loaded: {len(gene_cols)}")

    for c in gene_cols:
        expr[c] = pd.to_numeric(expr[c], errors="coerce")

    # ---------------------------
    # Patient-level expression:
    # ---------------------------
    expr_pat = (
        expr.groupby("patient_id", dropna=False, sort=False)
        .apply(lambda sub: pick_representative_sample(sub, gene_cols), include_groups=False)
        .reset_index()
    )

    # Keep only patient_id + gene columns
    keep_cols = ["patient_id"] + gene_cols
    expr_pat = expr_pat[keep_cols].copy()

    surv = load_survival(data_dir / "TCGA_PanCancer_Survival.txt")

    surv_pat = (
        surv.groupby("patient_id", dropna=False)
        .agg(
            stage=("stage", pick_max_stage),
            time=("time", first_nonnull),
            event=("event", first_nonnull),
            pfi_time=("pfi_time", first_nonnull),
            pfi_event=("pfi_event", first_nonnull),
            gender=("gender", first_nonnull),
            race=("race", first_nonnull),
            vital_status=("vital_status", first_nonnull),
            age_at_initial_pathologic_diagnosis=("age_at_initial_pathologic_diagnosis", first_nonnull),
            cancer_type_abbreviation=("cancer type abbreviation", first_nonnull),
        )
        .reset_index()
    )

    df = expr_pat.merge(surv_pat, on="patient_id", how="inner")

    # Rename clinical fields
    df = df.rename(
        columns={
            "age_at_initial_pathologic_diagnosis": "age_at_diagnosis",
            "cancer_type_abbreviation": "cancer_type",
        }
    )

    # Create *_expr columns 
    for ens, name in GENES.items():
        if ens in df.columns:
            df[f"{name}_expr"] = pd.to_numeric(df[ens], errors="coerce")

    df["time"] = pd.to_numeric(df["time"], errors="coerce")
    df["event"] = pd.to_numeric(df["event"], errors="coerce")

    # ---------------------------
    # PFI cutoff columns:
    # ---------------------------
    df["pfi_time"] = pd.to_numeric(df["pfi_time"], errors="coerce")
    df["pfi_event"] = pd.to_numeric(df["pfi_event"], errors="coerce")

    df["pfi_time_5yr"] = df["pfi_time"].clip(upper=PFI_CUTOFF_MONTHS)
    df["pfi_event_5yr"] = np.where(df["pfi_time"] > PFI_CUTOFF_MONTHS, 0, df["pfi_event"])

    df = df.dropna(subset=["stage", "CD24_expr", "time", "event"])

    out_path = out_dir / "tcga_analysis_table.csv"
    df.to_csv(out_path, index=False)

    print(f"\n✓ Wrote analysis table: {out_path}")
    print(f"Patients: {df['patient_id'].nunique()}")
    print(f"Stage counts: {df['stage'].value_counts(dropna=False).to_dict()}")
    print(f"PFI available: {df['pfi_time'].notna().any()}")
    print(f"PFI 5yr columns written: pfi_time_5yr, pfi_event_5yr (cutoff={PFI_CUTOFF_MONTHS} months)")


if __name__ == "__main__":
    main()