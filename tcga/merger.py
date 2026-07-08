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
    individual_layers: Dict[str, pd.DataFrame],
    out_dir: "pathlib.Path",
    study_id: str,
    data_types: list[str],
) -> None:
    """Persist the merged table, each individual layer, and a small manifest.

    Parameters
    ----------
    merged_df
        The result of :func:`merge_layers`.
    individual_layers
        Mapping from data‑type name to its original *wide* DataFrame (clinical
        already combined).  These are written to ``{type}.csv``.
    out_dir
        Output directory – will be created if it does not exist.
    study_id
        Identifier of the cBioPortal study – used in the manifest.
    data_types
        List of data‑type identifiers that were fetched (e.g. ``["clinical",
        "rnaseq", "cna", "mutations_wide"]``).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    # Write each layer
    for name, df in individual_layers.items():
        out_path = out_dir / f"{name}.csv"
        df.to_csv(out_path, index=False)

    # Write merged table
    merged_path = out_dir / "merged.csv"
    merged_df.to_csv(merged_path, index=False)

    # Manifest JSON – simple record of what was produced
    manifest = {
        "study_id": study_id,
        "fetched_at": pd.Timestamp.utcnow().isoformat() + "Z",
        "data_types": data_types,
        "samples": len(merged_df),
        "columns": len(merged_df.columns),
    }
    import json
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

# End of module
