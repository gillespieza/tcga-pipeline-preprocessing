"""
tcga/merger.py — Join the selected data layers on a unified sample identifier.

The module provides a single public function ``merge_layers`` that receives a
dictionary of DataFrames keyed by data‑type (e.g. ``{"clinical": df_clinical,
"rnaseq": df_rnaseq, ...}``).  It left‑joins each additional layer onto the
clinical base, ensuring that the final matrix contains *all* samples for which
clinical metadata exist.
"""

from __future__ import annotations

from typing import Dict
import pathlib
import pandas as pd

# ---------------------------------------------------------------------------
# Helper utilities
# ---------------------------------------------------------------------------

def _standardise_sample_id(df: pd.DataFrame, column: str = "SAMPLE_ID") -> pd.DataFrame:
    """Ensure the sample identifier column is present and uniformly formatted.

    TCGA barcodes can appear in various forms (e.g. ``TCGA-XX-XXXX-01A``).  For
    merging we keep the *full* barcode, as it is unique across all data types.
    """
    if column not in df.columns:
        raise KeyError(f"Expected column '{column}' not found in DataFrame")
    df = df.copy()
    df[column] = df[column].astype(str).str.upper()
    return df

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def merge_layers(layers: Dict[str, pd.DataFrame]) -> pd.DataFrame:
    """Merge a set of data layers on the ``SAMPLE_ID`` column.

    Parameters
    ----------
    layers
        Mapping from a layer name (``clinical``, ``rnaseq``, ``cna``, ``mutations_wide``
        or ``rppa``) to a pandas DataFrame that already contains a ``SAMPLE_ID``
        column.

    Returns
    -------
    pandas.DataFrame
        The merged wide table.  Columns are prefixed with the layer name to avoid
        collisions (e.g. ``rnaseq_GENE1``).  Rows correspond to samples present in
        the *clinical* dataframe; any missing values from optional layers are left
        as ``NaN``.
    """
    if "clinical" not in layers:
        raise ValueError("Clinical layer is required for merging")

    # Base dataframe – clinical (already contains SAMPLE_ID after combine_clinical)
    merged = _standardise_sample_id(layers["clinical"], "SAMPLE_ID")
    merged = merged.set_index("SAMPLE_ID")

    # Process each optional layer in a deterministic order
    for name in [k for k in layers.keys() if k != "clinical"]:
        df = _standardise_sample_id(layers[name], "SAMPLE_ID")
        df = df.set_index("SAMPLE_ID")
        # Prefix columns to keep them distinct after the join
        df = df.add_prefix(f"{name}_")
        merged = merged.join(df, how="left")

    merged = merged.reset_index()
    return merged

# ---------------------------------------------------------------------------
# Convenience function for writing CSVs and manifest
# ---------------------------------------------------------------------------

def write_outputs(
    merged_df: pd.DataFrame,
    raw_layers: Dict[str, pd.DataFrame],
    cleaned_layers: Dict[str, pd.DataFrame],
    out_dir: "pathlib.Path",
    study_id: str,
    data_types: list[str],
) -> None:
    """Persist outputs to two subdirectories:

    ``raw/``
        One CSV per fetched data type, exactly as returned by the API
        (clinical.csv, rnaseq.csv, cna.csv, rppa.csv, mutations_long.csv,
        mutations_wide.csv).

    ``processed/``
        clinical_cleaned.csv  — deduplicated, binary survival events, invalid
                                OS records removed.
        merged.csv            — all cleaned layers joined on SAMPLE_ID.
        manifest.json         — study ID, fetch timestamp, row/column counts.

    Parameters
    ----------
    merged_df
        The result of :func:`merge_layers` (built from cleaned_layers).
    raw_layers
        Mapping of data-type name -> raw DataFrame (as fetched from API).
    cleaned_layers
        Mapping of data-type name -> cleaned DataFrame.  The ``clinical``
        key holds the cleaned clinical table written as clinical_cleaned.csv.
    out_dir
        Study-level output directory – subdirs are created automatically.
    study_id
        cBioPortal study identifier – used in the manifest.
    data_types
        List of data-type identifiers that were fetched.
    """
    import json

    raw_dir       = out_dir / "raw"
    processed_dir = out_dir / "processed"
    raw_dir.mkdir(parents=True, exist_ok=True)
    processed_dir.mkdir(parents=True, exist_ok=True)

    # --- raw/ : one file per fetched layer, unaltered ---
    for name, df in raw_layers.items():
        df.to_csv(raw_dir / f"{name}.csv", index=False)

    # --- processed/ : cleaned layers + merged table ---
    for name, df in cleaned_layers.items():
        df.to_csv(processed_dir / f"{name}_cleaned.csv", index=False)

    merged_df.to_csv(processed_dir / "merged.csv", index=False)

    # Manifest lives in processed/ alongside the analysis-ready files
    manifest = {
        "study_id":   study_id,
        "fetched_at": pd.Timestamp.utcnow().isoformat() + "Z",
        "data_types": data_types,
        "raw_samples":       len(list(raw_layers.values())[0]) if raw_layers else 0,
        "processed_samples": len(merged_df),
        "processed_columns": len(merged_df.columns),
    }
    (processed_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )

# End of module
