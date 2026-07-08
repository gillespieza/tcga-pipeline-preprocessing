"""
tcga/loaders.py — Transform raw cBioPortal API responses into tidy DataFrames.

Each loader receives raw JSON records from the API (via tcga.api) and returns
a clean pandas DataFrame ready for merging.

Public functions:
  - build_clinical_df(sample_records, patient_records) -> DataFrame
  - build_molecular_df(records, value_col)             -> DataFrame (wide: samples × genes)
  - build_mutations_long(records)                      -> DataFrame (one row per variant)
  - build_mutations_wide(long_df)                      -> DataFrame (binary samples × genes)
"""

from __future__ import annotations

from typing import List, Any, Optional

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Clinical data helpers
# ---------------------------------------------------------------------------

def standardise_sample_id(sample_id: Any) -> str:
    """
    Standardise TCGA sample barcodes to a uniform 15-character hyphenated format
    (e.g., converts 'TCGA.KL.8323.01A' or 'TCGA-KL-8323-01A-11D' to 'TCGA-KL-8323-01')
    """
    if pd.isna(sample_id) or sample_id is None:
        return ""
    s = str(sample_id).strip().replace(".", "-").upper()
    if s.startswith("TCGA-") and len(s) > 15:
        return s[:15]
    return s


def parse_survival_status(val: Any) -> Optional[float]:
    """
    Map various clinical survival status string representations from cBioPortal
    to a binary 0.0 (censored/alive) or 1.0 (deceased/recurred/progressed).
    """
    if pd.isna(val) or val is None:
        return None
    val_str = str(val).strip().upper()
    # Matches codes starting with '1:' or exactly '1'
    if val_str.startswith("1:") or val_str == "1":
        return 1.0
    # Matches codes starting with '0:' or exactly '0'
    if val_str.startswith("0:") or val_str == "0":
        return 0.0
    # Fallback to keyword matching
    if any(word in val_str for word in ["DECEASED", "DEAD WITH TUMOR", "RECURRED", "PROGRESSION"]):
        return 1.0
    if any(word in val_str for word in ["LIVING", "ALIVE OR DEAD TUMOR FREE", "DISEASEFREE", "CENSORED"]):
        return 0.0
    return None


# ---------------------------------------------------------------------------
# Clinical data
# ---------------------------------------------------------------------------

def build_clinical_df(
    sample_records: List[dict],
    patient_records: List[dict],
) -> pd.DataFrame:
    """
    Pivot long-format clinical API records into a wide table with one row per
    sample, enriched with patient-level attributes. Returns the raw merged
    table without any survival cleaning or deduplication applied.

    API records look like:
        {"sampleId": "TCGA-XX-01", "clinicalAttributeId": "AGE", "value": "65"}
    """
    # --- Sample-level clinical data (pivot long -> wide) ---
    if not sample_records:
        raise ValueError("No sample-level clinical data returned from API")

    sample_df = pd.DataFrame(sample_records)
    sample_wide = sample_df.pivot_table(
        index="sampleId",
        columns="clinicalAttributeId",
        values="value",
        aggfunc="first",
    ).reset_index()
    sample_wide.columns.name = None
    sample_wide = sample_wide.rename(columns={"sampleId": "SAMPLE_ID"})

    # Extract PATIENT_ID from sample records if present
    if "patientId" in sample_df.columns:
        patient_map = (
            sample_df[["sampleId", "patientId"]]
            .drop_duplicates()
            .rename(columns={"sampleId": "SAMPLE_ID", "patientId": "PATIENT_ID"})
        )
        sample_wide = sample_wide.merge(patient_map, on="SAMPLE_ID", how="left")

    # --- Patient-level clinical data (pivot long -> wide) ---
    if patient_records:
        patient_df = pd.DataFrame(patient_records)
        patient_wide = patient_df.pivot_table(
            index="patientId",
            columns="clinicalAttributeId",
            values="value",
            aggfunc="first",
        ).reset_index()
        patient_wide.columns.name = None
        patient_wide = patient_wide.rename(columns={"patientId": "PATIENT_ID"})

        # Merge patient attributes onto sample table
        if "PATIENT_ID" in sample_wide.columns:
            # Avoid column collisions by suffixing duplicates
            sample_wide = sample_wide.merge(
                patient_wide,
                on="PATIENT_ID",
                how="left",
                suffixes=("", "_PATIENT"),
            )

    return sample_wide


def clean_clinical_df(raw_df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply deduplication and survival data cleaning to a raw clinical DataFrame
    produced by build_clinical_df().

    Steps:
      1. Standardise SAMPLE_ID.
      2. Deduplicate by SAMPLE_ID, then by PATIENT_ID (keep first sample per
         patient to ensure statistical independence in downstream models).
      3. Parse OS/DSS/DFS/PFS status strings into binary 0.0/1.0 event variables.
      4. Drop rows where Overall Survival (OS) is missing or has months <= 0.
      5. Set secondary endpoint (DSS/DFS/PFS) invalid/missing entries to NaN.

    Returns a new DataFrame; the input is not modified.
    """
    df = raw_df.copy()

    # Standardise sample IDs
    df["SAMPLE_ID"] = df["SAMPLE_ID"].apply(standardise_sample_id)
    df = df[df["SAMPLE_ID"] != ""]

    # --- 1. Deduplicate ---
    df = df.drop_duplicates(subset=["SAMPLE_ID"])

    if "PATIENT_ID" in df.columns:
        before_p_dedup = len(df)
        df = df.drop_duplicates(subset=["PATIENT_ID"])
        after_p_dedup = len(df)
        if before_p_dedup != after_p_dedup:
            print(f"\n  [INFO] Deduplicated patients: dropped {before_p_dedup - after_p_dedup} sample rows to ensure 1 sample per patient.")

    # --- 2. Clean Survival Records and Convert to Binary Event Variables ---
    # (status_col, months_col, is_primary)
    endpoints = [
        ("OS_STATUS",  "OS_MONTHS",  True),
        ("DSS_STATUS", "DSS_MONTHS", False),
        ("DFS_STATUS", "DFS_MONTHS", False),
        ("PFS_STATUS", "PFS_MONTHS", False),
    ]

    for status_col, months_col, is_primary in endpoints:
        if status_col in df.columns and months_col in df.columns:
            df[status_col] = df[status_col].apply(parse_survival_status)
            df[months_col] = pd.to_numeric(df[months_col], errors="coerce")

            invalid_mask = (
                df[months_col].isna() |
                (df[months_col] <= 0) |
                df[status_col].isna()
            )

            if is_primary:
                before_filter = len(df)
                df = df[~invalid_mask]
                after_filter = len(df)
                if before_filter != after_filter:
                    print(f"\n  [INFO] Removed {before_filter - after_filter} record(s) with missing or invalid overall survival (OS) data.")
            else:
                df.loc[invalid_mask, months_col] = np.nan
                df.loc[invalid_mask, status_col] = np.nan

    return df



# ---------------------------------------------------------------------------
# Molecular data (RNA-seq, CNA, RPPA)
# ---------------------------------------------------------------------------

def build_molecular_df(records: List[dict]) -> pd.DataFrame:
    """
    Pivot molecular data API records into a wide samples × genes matrix.

    API records look like:
        {"sampleId": "TCGA-XX-01", "entrezGeneId": 7157,
         "gene": {"hugoGeneSymbol": "TP53"}, "value": 1234.5}
    """
    if not records:
        return pd.DataFrame(columns=["SAMPLE_ID"])

    rows = []
    for r in records:
        gene_symbol = r.get("gene", {}).get("hugoGeneSymbol") or str(r.get("entrezGeneId", ""))
        rows.append({
            "SAMPLE_ID": r["sampleId"],
            "gene": gene_symbol,
            "value": r.get("value"),
        })

    long_df = pd.DataFrame(rows)
    long_df["value"] = pd.to_numeric(long_df["value"], errors="coerce")

    # Pivot to wide: rows = samples, columns = genes
    wide_df = long_df.pivot_table(
        index="SAMPLE_ID",
        columns="gene",
        values="value",
        aggfunc="first",
    ).reset_index()
    wide_df.columns.name = None

    return wide_df


def clean_rnaseq_df(rnaseq_df: pd.DataFrame) -> pd.DataFrame:
    """
    Deduplicate RNA-seq sample records (ensuring one sample per patient)
    and remove columns (genes) with missing values.
    """
    if rnaseq_df.empty or "SAMPLE_ID" not in rnaseq_df.columns:
        return rnaseq_df

    df = rnaseq_df.copy()

    # Standardise sample IDs
    df["SAMPLE_ID"] = df["SAMPLE_ID"].apply(standardise_sample_id)
    df = df[df["SAMPLE_ID"] != ""]

    # --- 1. Deduplicate by Patient ID ---
    # Extract PATIENT_ID from SAMPLE_ID (first 12 chars for TCGA barcode)
    df["PATIENT_ID"] = df["SAMPLE_ID"].apply(
        lambda x: "-".join(x.split("-")[:3]) if isinstance(x, str) and x.startswith("TCGA-") else x
    )

    before_p_dedup = len(df)
    df = df.drop_duplicates(subset=["PATIENT_ID"])
    after_p_dedup = len(df)
    if before_p_dedup != after_p_dedup:
        print(f"\n  [INFO] RNA-seq deduplicated patients: dropped {before_p_dedup - after_p_dedup} sample rows to ensure 1 sample per patient.")

    df = df.drop(columns=["PATIENT_ID"])

    # --- 2. Remove Missing Values (Genes/Columns with NaN) ---
    gene_cols = [c for c in df.columns if c != "SAMPLE_ID"]
    missing_counts = df[gene_cols].isna().sum()
    cols_with_nans = missing_counts[missing_counts > 0].index.tolist()
    if cols_with_nans:
        df = df.drop(columns=cols_with_nans)
        print(f"\n  [INFO] RNA-seq: Removed {len(cols_with_nans)} gene(s) with missing values.")

    return df


def clean_cna_df(cna_df: pd.DataFrame) -> pd.DataFrame:
    """
    Deduplicate CNA sample records (ensuring one sample per patient)
    and remove columns (genes) with missing values.
    """
    if cna_df.empty or "SAMPLE_ID" not in cna_df.columns:
        return cna_df

    df = cna_df.copy()

    # Standardise sample IDs
    df["SAMPLE_ID"] = df["SAMPLE_ID"].apply(standardise_sample_id)
    df = df[df["SAMPLE_ID"] != ""]

    # --- 1. Deduplicate by Patient ID ---
    # Extract PATIENT_ID from SAMPLE_ID (first 12 chars for TCGA barcode)
    df["PATIENT_ID"] = df["SAMPLE_ID"].apply(
        lambda x: "-".join(x.split("-")[:3]) if isinstance(x, str) and x.startswith("TCGA-") else x
    )

    before_p_dedup = len(df)
    df = df.drop_duplicates(subset=["PATIENT_ID"])
    after_p_dedup = len(df)
    if before_p_dedup != after_p_dedup:
        print(f"\n  [INFO] CNA deduplicated patients: dropped {before_p_dedup - after_p_dedup} sample rows to ensure 1 sample per patient.")

    df = df.drop(columns=["PATIENT_ID"])

    # --- 2. Remove Missing Values (Genes/Columns with NaN) ---
    gene_cols = [c for c in df.columns if c != "SAMPLE_ID"]
    missing_counts = df[gene_cols].isna().sum()
    cols_with_nans = missing_counts[missing_counts > 0].index.tolist()
    if cols_with_nans:
        df = df.drop(columns=cols_with_nans)
        print(f"\n  [INFO] CNA: Removed {len(cols_with_nans)} gene(s) with missing values.")

    return df


# ---------------------------------------------------------------------------
# Mutation data
# ---------------------------------------------------------------------------

def build_mutations_long(records: List[dict]) -> pd.DataFrame:
    """
    Build a long-format mutations table (one row per variant) from API records.

    Keeps the most commonly used MAF-like columns.
    """
    if not records:
        return pd.DataFrame(columns=["SAMPLE_ID", "HUGO_SYMBOL"])

    rows = []
    for r in records:
        rows.append({
            "SAMPLE_ID":              r.get("sampleId", ""),
            "HUGO_SYMBOL":            r.get("gene", {}).get("hugoGeneSymbol", ""),
            "ENTREZ_GENE_ID":         r.get("entrezGeneId", ""),
            "VARIANT_CLASSIFICATION": r.get("mutationType", ""),
            "VARIANT_TYPE":           r.get("variantType", ""),
            "MUTATION_STATUS":        r.get("mutationStatus", ""),
            "CHR":                    r.get("chr", ""),
            "START_POSITION":         r.get("startPosition", ""),
            "END_POSITION":           r.get("endPosition", ""),
            "REF_ALLELE":             r.get("referenceAllele", ""),
            "ALT_ALLELE":             r.get("variantAllele", ""),
            "PROTEIN_CHANGE":         r.get("proteinChange", ""),
            "KEYWORD":                r.get("keyword", ""),
            "NCBI_BUILD":             r.get("ncbiBuild", ""),
        })

    return pd.DataFrame(rows)


def build_mutations_wide(long_df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a binary sample × gene mutation matrix from the long-format table.

    A value of 1 means the gene has at least one mutation in that sample.
    """
    if long_df.empty or "SAMPLE_ID" not in long_df.columns:
        return pd.DataFrame(columns=["SAMPLE_ID"])

    binary = (
        long_df[["SAMPLE_ID", "HUGO_SYMBOL"]]
        .drop_duplicates()
        .assign(MUT=1)
        .pivot_table(index="SAMPLE_ID", columns="HUGO_SYMBOL", values="MUT", aggfunc="first")
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    binary.columns.name = None
    return binary


def clean_mutations_df(
    long_df: pd.DataFrame,
    print_fn: Optional[Any] = None,
) -> pd.DataFrame:
    """
    Clean somatic mutation records (long format):
      1. Standardise SAMPLE_ID.
      2. Filter out missing/empty IDs.
      3. Filter for non-silent/functional mutations.
      4. Deduplicate patients to ensure 1 sample per patient.

    Args:
        long_df:  Raw long-format mutation DataFrame.
        print_fn: Optional callable for progress output (e.g. rich's rprint).
                  When provided, a status line is emitted after each step.
    """
    _log = print_fn if callable(print_fn) else lambda *a, **kw: None

    if long_df.empty:
        return long_df

    df = long_df.copy()
    n_start = len(df)

    # 1. Standardise sample IDs
    _log(f"    [dim]Standardising sample IDs ({n_start:,} variants)...[/]")
    df["SAMPLE_ID"] = df["SAMPLE_ID"].apply(standardise_sample_id)

    # 2. Filter out missing/empty IDs or HUGO_SYMBOL
    before_filter = len(df)
    df = df[df["SAMPLE_ID"].str.strip() != ""]
    df = df[df["HUGO_SYMBOL"].notna() & (df["HUGO_SYMBOL"].str.strip() != "")]
    dropped_missing = before_filter - len(df)
    if dropped_missing:
        _log(f"    [dim]Dropped {dropped_missing:,} record(s) with missing ID or gene symbol.[/]")

    # 3. Filter for non-silent/functional mutations
    non_silent = {
        "FRAME_SHIFT_DEL", "FRAME_SHIFT_INS", "IN_FRAME_DEL", "IN_FRAME_INS",
        "MISSENSE_MUTATION", "NONSENSE_MUTATION", "SPLICE_SITE",
        "TRANSLATION_START_SITE", "NONSTOP_MUTATION"
    }
    before_ns = len(df)
    df = df[df["VARIANT_CLASSIFICATION"].astype(str).str.upper().isin(non_silent)]
    dropped_silent = before_ns - len(df)
    _log(
        f"    [dim]Filtered silent/non-functional variants: "
        f"{before_ns:,} → {len(df):,} "
        f"(removed {dropped_silent:,} silent/intronic/UTR).[/]"
    )

    # 4. Deduplicate patients (keep 1 sample per patient)
    df["PATIENT_ID"] = df["SAMPLE_ID"].apply(
        lambda x: "-".join(x.split("-")[:3]) if isinstance(x, str) and x.startswith("TCGA-") else x
    )
    n_patients_before = df["PATIENT_ID"].nunique()
    n_samples_before  = df["SAMPLE_ID"].nunique()

    # Keep first unique sample_id per patient
    patient_sample_map = (
        df[["PATIENT_ID", "SAMPLE_ID"]]
        .drop_duplicates()
        .groupby("PATIENT_ID")
        .first()
        .reset_index()
    )
    df = df[df["SAMPLE_ID"].isin(patient_sample_map["SAMPLE_ID"])]

    n_samples_after = df["SAMPLE_ID"].nunique()
    dropped_dup_samples = n_samples_before - n_samples_after
    if dropped_dup_samples:
        _log(
            f"    [dim]Patient deduplication: {n_samples_before} samples → "
            f"{n_samples_after} ({dropped_dup_samples} duplicate sample(s) removed, "
            f"{n_patients_before} unique patients retained).[/]"
        )
    else:
        _log(f"    [dim]Patient deduplication: {n_samples_after} samples, no duplicates found.[/]")

    # Drop temporary column
    df = df.drop(columns=["PATIENT_ID"])

    return df
