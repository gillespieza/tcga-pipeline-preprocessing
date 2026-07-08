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

from typing import List

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Clinical data
# ---------------------------------------------------------------------------

def build_clinical_df(
    sample_records: List[dict],
    patient_records: List[dict],
) -> pd.DataFrame:
    """
    Pivot long-format clinical API records into a wide table with one row per
    sample, enriched with patient-level attributes.

    API records look like:
        {"sampleId": "TCGA-XX-01", "clinicalAttributeId": "AGE", "value": "65"}
    """
    # --- Sample-level clinical data (pivot long → wide) ---
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

    # --- Patient-level clinical data (pivot long → wide) ---
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
