#!/usr/bin/env python3
"""
TCGA colorectal cancer survival analysis
- Stratified by CD24 and Siglec-10 expression

Inputs (from UCSC Xena - check README file):
- Gene expression matrices
    - TCGA-COAD.star_fpkm-uq.tsv.gz
    - TCGA-READ.star_fpkm-uq.tsv.gz
- Clinical data
    - TCGA-COAD.clinical.tsv(.gz)
    - TCGA-READ.clinical.tsv(.gz)
- Survival data
    - TCGA_PanCancer_Survival.txt

Outputs:
- TCGA_Patient_Data.csv
- TCGA_CD24_Survival_Summary.csv
- TCGA_CD24_Siglec10_Combined_Survival_Summary.csv
- Kaplan-Meier Plots (OS + PFI)

"""

from __future__ import annotations
import argparse
import gzip
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test


# ----------------------------
# Helpers
# ----------------------------
def read_table(path: Path, sep: str = "\t") -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing input file: {path}")
    if path.suffix == ".gz":
        return pd.read_csv(path, sep=sep, compression="gzip", low_memory=False)
    return pd.read_csv(path, sep=sep, low_memory=False)


def read_xena_expression(path: Path) -> pd.DataFrame:

    df = read_table(path, sep="\t")
    first_col = df.columns[0]
    df = df.rename(columns={first_col: "gene"})
    df["gene"] = df["gene"].astype(str)
    df = df.set_index("gene")
    return df


def tcga_barcode_to_patient(barcode: str) -> str:
    return str(barcode)[:12]


def normalize_stage(raw: str) -> Optional[str]:
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return None
    s = str(raw).strip()
    if s == "" or s.lower() in {"na", "nan", "not reported"}:
        return None

    m = re.search(r"stage\s*([iv]+)", s.lower())
    if not m:
        return None
    roman = m.group(1).upper()
    if roman in {"I", "II", "III", "IV"}:
        return roman
    return None


def safe_log2_fpkm_uq(x: pd.Series) -> pd.Series:
    x = pd.to_numeric(x, errors="coerce")
    x = x.fillna(0.0)
    x = x.clip(lower=0.0)
    return np.log2(x + 1.0)


def median_group(series: pd.Series, high_label="High", low_label="Low") -> Tuple[pd.Series, float]:
    s = pd.to_numeric(series, errors="coerce")
    med = float(np.nanmedian(s.values))
    grp = np.where(s >= med, high_label, low_label)
    return pd.Series(grp, index=series.index), med


def km_plot_2group(df: pd.DataFrame, time_col: str, event_col: str, group_col: str,
                   title: str, outpath: Optional[Path] = None) -> None:
    plt.figure(figsize=(5, 4))
    kmf = KaplanMeierFitter()

    for g in ["Low", "High"]:
        sub = df[df[group_col] == g]
        if sub.empty:
            continue
        kmf.fit(sub[time_col], event_observed=sub[event_col], label=g)
        kmf.plot(ci_show=False)

    plt.title(title)
    plt.xlabel("Time (days)")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    if outpath is not None:
        outpath.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close()


def logrank_2group(df: pd.DataFrame, time_col: str, event_col: str, group_col: str) -> float:
    a = df[df[group_col] == "Low"]
    b = df[df[group_col] == "High"]
    if a.empty or b.empty:
        return np.nan
    res = logrank_test(a[time_col], b[time_col], event_observed_A=a[event_col], event_observed_B=b[event_col])
    return float(res.p_value)


def logrank_4group(df: pd.DataFrame, time_col: str, event_col: str, group_col: str) -> float:
    sub = df.dropna(subset=[time_col, event_col, group_col]).copy()
    if sub[group_col].nunique() < 2:
        return np.nan
    res = multivariate_logrank_test(sub[time_col], sub[group_col], sub[event_col])
    return float(res.p_value)


# ----------------------------
# Main
# ----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_dir", type=Path, required=True,
                    help="Folder containing the downloaded Xena files.")
    ap.add_argument("--out_dir", type=Path, required=True,
                    help="Folder to write outputs (CSVs and optional plots).")
    ap.add_argument("--coad_expr", type=str, default="TCGA-COAD.star_fpkm-uq.tsv.gz")
    ap.add_argument("--read_expr", type=str, default="TCGA-READ.star_fpkm-uq.tsv.gz")
    ap.add_argument("--coad_clin", type=str, default="TCGA-COAD.clinical.tsv")
    ap.add_argument("--read_clin", type=str, default="TCGA-READ.clinical.tsv")
    ap.add_argument("--survival", type=str, default="TCGA_PanCancer_Survival.txt")
    ap.add_argument("--out_patient", type=str, default="TCGA_Patient_Data.csv")
    ap.add_argument("--out_cd24", type=str, default="TCGA_CD24_Survival_Summary.csv")
    ap.add_argument("--out_combo", type=str, default="TCGA_CD24_Siglec10_Combined_Survival_Summary.csv")

    # TAM gene set
    ap.add_argument(
        "--tam_genes",
        type=str,
        default="SIGLEC10,C1QA,C1QB,C1QC,APOE,LST1,CSF1R,TYROBP,FCER1G,AIF1,CD68,MSR1",
        help="Comma-separated TAM gene list used to compute Siglec10_TAM_Score.",
    )
    ap.add_argument("--make_plots", action="store_true", help="Also export KM plots (optional).")

    args = ap.parse_args()

    in_dir = args.input_dir
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- load expression ----
    expr_coad = read_xena_expression(in_dir / args.coad_expr)
    expr_read = read_xena_expression(in_dir / args.read_expr)

    # concatenate cohorts by columns (samples)
    expr = pd.concat([expr_coad, expr_read], axis=1)

    # ---- load clinical ----
    clin_coad = read_table(in_dir / args.coad_clin, sep="\t")
    clin_read = read_table(in_dir / args.read_clin, sep="\t")
    clin = pd.concat([clin_coad, clin_read], ignore_index=True)

    # Xena clinical usually includes "sample" column
    if "sample" not in clin.columns:
        # try common alternatives
        for cand in ["Sample", "sampleID", "submitter_id"]:
            if cand in clin.columns:
                clin = clin.rename(columns={cand: "sample"})
                break
    if "sample" not in clin.columns:
        raise ValueError(f"Clinical file missing sample id column. Columns: {list(clin.columns)}")

    clin["patient_id"] = clin["sample"].map(tcga_barcode_to_patient)

    # stage column
    stage_col = "ajcc_pathologic_stage.diagnoses"
    if stage_col not in clin.columns:
        # fallback search
        candidates = [c for c in clin.columns if "pathologic_stage" in c.lower() or "ajcc" in c.lower() and "stage" in c.lower()]
        if not candidates:
            raise ValueError(f"Could not find a stage column in clinical. Columns: {list(clin.columns)}")
        stage_col = candidates[0]

    clin["Stage"] = clin[stage_col].map(normalize_stage)

    # ---- survival ----
    surv = read_table(in_dir / args.survival, sep="\t")
    # survival table has sample IDs; harmonize to patient-level
    if "sample" not in surv.columns:
        # first col is "sample"
        surv = surv.rename(columns={surv.columns[0]: "sample"})
    surv["patient_id"] = surv["sample"].map(tcga_barcode_to_patient)

    # Keep patient-level unique entries (first occurrence)
    surv = surv.sort_values("sample").drop_duplicates("patient_id", keep="first")

    # ---- Build patient table ----
    # Map cohort by sample barcode: COAD vs READ
    cohort_map = {}
    for s in expr_coad.columns:
        cohort_map[tcga_barcode_to_patient(s)] = "COAD"
    for s in expr_read.columns:
        cohort_map[tcga_barcode_to_patient(s)] = "READ"

    # Extract CD24 and SIGLEC10 vectors
    for g in ["CD24", "SIGLEC10"]:
        if g not in expr.index:
            raise ValueError(
                f"Gene '{g}' not found in expression matrix rownames. "
                f"Tip: check expr.index[:20] to see whether genes are symbols or Ensembl IDs."
            )

    cd24 = safe_log2_fpkm_uq(expr.loc["CD24"])
    sig10 = safe_log2_fpkm_uq(expr.loc["SIGLEC10"])

    # TAM score
    tam_genes = [x.strip() for x in args.tam_genes.split(",") if x.strip()]
    missing_tam = [g for g in tam_genes if g not in expr.index]
    if missing_tam:
        raise ValueError(
            "These TAM genes are missing from expression matrix rownames:\n"
            + ", ".join(missing_tam)
            + "\n\nEdit --tam_genes or add a symbol/Ensembl mapping step."
        )

    tam_mat = expr.loc[tam_genes].apply(safe_log2_fpkm_uq, axis=1)
    tam_score = tam_mat.mean(axis=0)

    # Patient IDs per sample column
    patient_ids = pd.Index([tcga_barcode_to_patient(x) for x in expr.columns], name="patient_id")

    df_pat = pd.DataFrame(
        {
            "patient_id": patient_ids,
            "cohort": [cohort_map.get(pid, "COAD/READ") for pid in patient_ids],
            "CD24_Expression_Log2": cd24.values,
            "SIGLEC10_Expression_Log2": sig10.values,
            "Siglec10_TAM_Score": tam_score.values,
        }
    ).drop_duplicates("patient_id", keep="first")

    # add stage (patient-level)
    clin_pat = clin.sort_values("sample").drop_duplicates("patient_id", keep="first")[["patient_id", "Stage"]]
    df_pat = df_pat.merge(clin_pat, on="patient_id", how="left")

    # add survival endpoints (patient-level)
    keep_surv = ["patient_id", "OS", "OS.time", "PFI", "PFI.time"]
    missing_cols = [c for c in keep_surv if c not in surv.columns]
    if missing_cols:
        raise ValueError(f"Survival table missing columns: {missing_cols}. Columns: {list(surv.columns)}")

    df_pat = df_pat.merge(surv[keep_surv], on="patient_id", how="left")

    # group splits
    df_pat["CD24_Group"], cd24_med = median_group(df_pat["CD24_Expression_Log2"], high_label="High", low_label="Low")
    df_pat["Siglec10_TAM_Group"], tam_med = median_group(df_pat["Siglec10_TAM_Score"], high_label="High", low_label="Low")

    def combined_group(row) -> str:
        return f"{row['CD24_Group']}/{row['Siglec10_TAM_Group']}"

    df_pat["Combined_Group"] = df_pat.apply(combined_group, axis=1)

    # Save patient table (this is your “main output”)
    df_pat = df_pat.sort_values(["Stage", "patient_id"]).reset_index(drop=True)
    df_pat.to_csv(out_dir / args.out_patient, index=False)

    # ---- Survival summaries ----
    def summarize_cd24(df_stage: pd.DataFrame, endpoint: str) -> Dict:
        time_col = f"{endpoint}.time"
        event_col = endpoint
        sub = df_stage.dropna(subset=[time_col, event_col, "CD24_Group"]).copy()
        sub[time_col] = pd.to_numeric(sub[time_col], errors="coerce")
        sub[event_col] = pd.to_numeric(sub[event_col], errors="coerce")
        sub = sub.dropna(subset=[time_col, event_col])

        p = logrank_2group(sub, time_col, event_col, "CD24_Group")

        out = {
            "Endpoint": endpoint,
            "Stage": sub["Stage"].iloc[0] if len(sub) else df_stage["Stage"].iloc[0],
            "N (tumor samples)": int(len(sub)),
            "Events": int(sub[event_col].sum()) if len(sub) else 0,
            "Log-rank p": p,
            "CD24_Low": int((sub["CD24_Group"] == "Low").sum()),
            "CD24_High": int((sub["CD24_Group"] == "High").sum()),
        }
        return out

    def summarize_combo(df_stage: pd.DataFrame, endpoint: str) -> Dict:
        time_col = f"{endpoint}.time"
        event_col = endpoint
        sub = df_stage.dropna(subset=[time_col, event_col, "Combined_Group"]).copy()
        sub[time_col] = pd.to_numeric(sub[time_col], errors="coerce")
        sub[event_col] = pd.to_numeric(sub[event_col], errors="coerce")
        sub = sub.dropna(subset=[time_col, event_col])

        p4 = logrank_4group(sub, time_col, event_col, "Combined_Group")

        # counts per combined group
        counts = sub["Combined_Group"].value_counts().to_dict()
        out = {
            "Endpoint": endpoint,
            "Stage": sub["Stage"].iloc[0] if len(sub) else df_stage["Stage"].iloc[0],
            "N (tumor samples)": int(len(sub)),
            "Events": int(sub[event_col].sum()) if len(sub) else 0,
            "Log-rank p (4-group)": p4,
            "Low/Low": int(counts.get("Low/Low", 0)),
            "Low/High": int(counts.get("Low/High", 0)),
            "High/Low": int(counts.get("High/Low", 0)),
            "High/High": int(counts.get("High/High", 0)),
        }
        return out

    # stage-wise summaries (I–IV)
    stages = ["I", "II", "III", "IV"]
    cd24_rows = []
    combo_rows = []

    for st in stages:
        df_st = df_pat[df_pat["Stage"] == st].copy()
        if df_st.empty:
            continue
        for ep in ["OS", "PFI"]:
            cd24_rows.append(summarize_cd24(df_st, ep))
            combo_rows.append(summarize_combo(df_st, ep))

        if args.make_plots:
            # optional plots (not required for your repo)
            for ep in ["OS", "PFI"]:
                km_plot_2group(
                    df_st, time_col=f"{ep}.time", event_col=ep, group_col="CD24_Group",
                    title=f"{ep} | Stage {st} | CD24 High vs Low",
                    outpath=out_dir / "plots" / f"KM_{ep}_Stage{st}_CD24.png",
                )

    df_cd24 = pd.DataFrame(cd24_rows).sort_values(["Endpoint", "Stage"]).reset_index(drop=True)
    df_combo = pd.DataFrame(combo_rows).sort_values(["Endpoint", "Stage"]).reset_index(drop=True)

    df_cd24.to_csv(out_dir / args.out_cd24, index=False)
    df_combo.to_csv(out_dir / args.out_combo, index=False)

    print("Saved:")
    print(" -", out_dir / args.out_patient)
    print(" -", out_dir / args.out_cd24)
    print(" -", out_dir / args.out_combo)
    print("\nMedians used:")
    print(f" - CD24 median (log2): {cd24_med:.6f}")
    print(f" - TAM median  (log2): {tam_med:.6f}")


if __name__ == "__main__":
    main()