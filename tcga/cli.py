"""
tcga/cli.py — Interactive and scriptable CLI for the TCGA preprocessing pipeline.

Uses the cBioPortal REST API directly to fetch all data (no tarball download).
Supports:
  1. Interactive — fuzzy-search for a study, checkbox data-type selection (questionary).
  2. CLI args   — fully scriptable via argparse flags.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import questionary
from rich import print as rprint
from rich.panel import Panel
from rich.table import Table

from .api import (
    detect_available_data_types,
    get_clinical_data,
    get_molecular_data,
    get_mutations,
    get_sample_ids,
    get_study,
    list_studies,
    resolve_profile_id,
)
from .loaders import (
    build_clinical_df,
    build_molecular_df,
    build_mutations_long,
    build_mutations_wide,
)
from .merger import merge_layers, write_outputs


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="tcga_fetch.py",
        description="Download and preprocess TCGA multi-omics data from cBioPortal.",
    )
    parser.add_argument(
        "--study",
        help="cBioPortal study ID (e.g. kirc_tcga_pan_can_atlas_2018). "
             "If omitted, launches interactive selection.",
    )
    parser.add_argument(
        "--data",
        nargs="+",
        choices=["clinical", "rnaseq", "cna", "mutations", "rppa"],
        help="Data types to fetch. Clinical is always included as the merge base.",
    )
    parser.add_argument(
        "--output",
        default="output",
        help="Root output directory (default: ./output).",
    )
    parser.add_argument(
        "--log2",
        action="store_true",
        help="Apply log2(x+1) transform to RNA-seq values.",
    )
    parser.add_argument(
        "--no-merge",
        action="store_true",
        help="Write individual layer CSVs only — skip merged.csv.",
    )
    parser.add_argument(
        "--list-studies",
        action="store_true",
        help="Print all available cBioPortal studies and exit.",
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Interactive prompts
# ---------------------------------------------------------------------------

def _list_studies_and_exit() -> None:
    """Print a rich table of all available studies, then exit."""
    rprint("[bold]Fetching study list from cBioPortal...[/]")
    studies = list_studies()

    table = Table(title="Available cBioPortal Studies", show_lines=False)
    table.add_column("Study ID", style="cyan", no_wrap=True)
    table.add_column("Name", style="magenta")
    table.add_column("Cancer Type", style="green")
    table.add_column("Samples", style="yellow", justify="right")

    for s in sorted(studies, key=lambda x: x["studyId"]):
        table.add_row(
            s["studyId"],
            s["name"],
            s["cancerType"],
            str(s["allSampleCount"]),
        )
    rprint(table)
    sys.exit(0)


def _interactive_study_selection() -> str:
    """Fuzzy-search study picker powered by questionary."""
    rprint("[bold]Fetching study list from cBioPortal...[/]")
    studies = list_studies()

    choices = [
        f"{s['studyId']}  —  {s['name']}  ({s['cancerType']})"
        for s in sorted(studies, key=lambda x: x["studyId"])
    ]

    chosen = questionary.autocomplete(
        "Search for a study (type to filter):",
        choices=choices,
        validate=lambda x: bool(x),
        ignore_case=True,
        match_middle=True,
    ).ask()

    if not chosen:
        sys.exit("No study selected — exiting.")

    # Extract the study ID (everything before the first whitespace)
    study_id = chosen.split()[0]
    return study_id


def _interactive_data_selection(study_id: str) -> List[str]:
    """Multi-select checkbox for data types, greying out unavailable ones."""
    rprint("[bold]Checking available data types...[/]")
    avail = detect_available_data_types(study_id)

    type_labels = {
        "clinical":  "Clinical (patient + sample)",
        "rnaseq":    "RNA-seq expression (mRNA)",
        "cna":       "Copy-number alterations (CNA)",
        "mutations": "Somatic mutations (MAF)",
        "rppa":      "RPPA proteomics",
    }

    choices = []
    for key, label in type_labels.items():
        is_available = avail.get(key, False)
        choices.append(
            questionary.Choice(
                title=label,
                value=key,
                checked=is_available,
                disabled="Not available for this study" if not is_available else None,
            )
        )

    selected = questionary.checkbox(
        "Select data types to fetch:",
        choices=choices,
    ).ask()

    if not selected:
        sys.exit("No data types selected — exiting.")

    return selected


def _interactive_log2_prompt() -> bool:
    """Ask whether to apply log2 transform to RNA-seq."""
    return questionary.confirm(
        "Apply log2(x+1) transform to RNA-seq values?",
        default=False,
    ).ask()


# ---------------------------------------------------------------------------
# Log2 transform helper
# ---------------------------------------------------------------------------

def _apply_log2_transform(df: pd.DataFrame) -> pd.DataFrame:
    """Apply log2(x+1) to all numeric columns except SAMPLE_ID."""
    df = df.copy()
    numeric_cols = df.select_dtypes(include=["number"]).columns
    for col in numeric_cols:
        df[col] = np.log2(df[col].clip(lower=0) + 1)
    return df


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> None:
    args = _parse_args(argv)

    # Banner
    rprint(Panel.fit(
        "[bold white]TCGA Data Preprocessing Pipeline[/]",
        border_style="cyan",
    ))

    # --list-studies: print and exit
    if args.list_studies:
        _list_studies_and_exit()

    # ── 1. Study selection ────────────────────────────────────────────
    study_id = args.study or _interactive_study_selection()

    # Validate the study ID exists
    meta = get_study(study_id)
    if meta is None:
        rprint(f"[bold red]Error:[/] Study '{study_id}' not found on cBioPortal.")
        sys.exit(1)
    rprint(f"[bold cyan]Study:[/] {study_id}  —  {meta.get('name', '')}")

    # ── 2. Data type selection ────────────────────────────────────────
    if args.data:
        data_requested = list(set(args.data) | {"clinical"})
    else:
        data_requested = _interactive_data_selection(study_id)
        if "clinical" not in data_requested:
            data_requested.insert(0, "clinical")

    rprint(f"[bold]Data types:[/] {', '.join(sorted(data_requested))}")

    # ── 3. Log2 decision ─────────────────────────────────────────────
    apply_log2 = args.log2
    if not args.study and "rnaseq" in data_requested and not args.log2:
        apply_log2 = _interactive_log2_prompt()

    # ── 4. Get sample IDs for the study ───────────────────────────────
    rprint()
    rprint("  Fetching sample list...", end=" ")
    sample_ids = get_sample_ids(study_id)
    rprint(f"[green]OK[/] {len(sample_ids)} samples")

    # ── 5. Fetch each requested layer via REST API ────────────────────
    layers: dict[str, pd.DataFrame] = {}
    individual_files: dict[str, pd.DataFrame] = {}

    # --- Clinical ---
    if "clinical" in data_requested:
        rprint("  Fetching clinical data...")
        try:
            sample_clinical = get_clinical_data(study_id, "SAMPLE")
            patient_clinical = get_clinical_data(study_id, "PATIENT")
            clinical_df = build_clinical_df(sample_clinical, patient_clinical)
            layers["clinical"] = clinical_df
            individual_files["clinical"] = clinical_df
            rprint(f"  [green]OK[/] {len(clinical_df)} samples, {len(clinical_df.columns)} attributes")
        except Exception as exc:
            rprint(f"  [red]ERR[/] {exc}")
            sys.exit(1)

    # --- RNA-seq ---
    if "rnaseq" in data_requested:
        rprint("  Fetching RNA-seq data...")
        profile_id = resolve_profile_id(study_id, "rnaseq")
        if profile_id:
            try:
                records = get_molecular_data(profile_id, sample_ids)
                rnaseq_df = build_molecular_df(records)
                if apply_log2:
                    rnaseq_df = _apply_log2_transform(rnaseq_df)
                    rprint(f"  [green]OK[/] {len(rnaseq_df)} samples x "
                           f"{len(rnaseq_df.columns) - 1} genes [log2(x+1)]")
                else:
                    rprint(f"  [green]OK[/] {len(rnaseq_df)} samples x "
                           f"{len(rnaseq_df.columns) - 1} genes [raw]")
                layers["rnaseq"] = rnaseq_df
                individual_files["rnaseq"] = rnaseq_df
            except Exception as exc:
                rprint(f"  [yellow]WARN[/] RNA-seq fetch failed: {exc}")
        else:
            rprint("  [yellow]WARN[/] No RNA-seq profile found — skipping")

    # --- CNA ---
    if "cna" in data_requested:
        rprint("  Fetching CNA data...")
        profile_id = resolve_profile_id(study_id, "cna")
        if profile_id:
            try:
                records = get_molecular_data(profile_id, sample_ids)
                cna_df = build_molecular_df(records)
                rprint(f"  [green]OK[/] {len(cna_df)} samples x "
                       f"{len(cna_df.columns) - 1} genes")
                layers["cna"] = cna_df
                individual_files["cna"] = cna_df
            except Exception as exc:
                rprint(f"  [yellow]WARN[/] CNA fetch failed: {exc}")
        else:
            rprint("  [yellow]WARN[/] No CNA profile found — skipping")

    # --- RPPA ---
    if "rppa" in data_requested:
        rprint("  Fetching RPPA data...")
        profile_id = resolve_profile_id(study_id, "rppa")
        if profile_id:
            try:
                records = get_molecular_data(profile_id, sample_ids)
                rppa_df = build_molecular_df(records)
                rprint(f"  [green]OK[/] {len(rppa_df)} samples x "
                       f"{len(rppa_df.columns) - 1} proteins")
                layers["rppa"] = rppa_df
                individual_files["rppa"] = rppa_df
            except Exception as exc:
                rprint(f"  [yellow]WARN[/] RPPA fetch failed: {exc}")
        else:
            rprint("  [yellow]WARN[/] No RPPA profile found — skipping")

    # --- Mutations ---
    if "mutations" in data_requested:
        rprint("  Fetching mutation data...")
        profile_id = resolve_profile_id(study_id, "mutations")
        if profile_id:
            try:
                records = get_mutations(profile_id, sample_ids)
                long_df = build_mutations_long(records)
                wide_df = build_mutations_wide(long_df)
                rprint(f"  [green]OK[/] {len(long_df)} variants -> "
                       f"{len(wide_df)} samples x {len(wide_df.columns) - 1} genes (wide)")
                individual_files["mutations_long"] = long_df
                individual_files["mutations_wide"] = wide_df
                layers["mutations_wide"] = wide_df
            except Exception as exc:
                rprint(f"  [yellow]WARN[/] Mutations fetch failed: {exc}")
        else:
            rprint("  [yellow]WARN[/] No mutation profile found — skipping")

    # ── 6. Merge and write outputs ────────────────────────────────────
    out_dir = Path(args.output) / study_id
    rprint()

    if not args.no_merge and "clinical" in layers:
        merge_eligible = {k: v for k, v in layers.items()}
        merged_df = merge_layers(merge_eligible)
        write_outputs(merged_df, individual_files, out_dir, study_id, sorted(data_requested))
        rprint(f"  [green]OK[/] Merged: {len(merged_df)} samples x {len(merged_df.columns)} columns")
    else:
        out_dir.mkdir(parents=True, exist_ok=True)
        for name, df in individual_files.items():
            filename = "clinical_cleaned.csv" if name == "clinical" else f"{name}.csv"
            df.to_csv(out_dir / filename, index=False)

    rprint()
    rprint(f"[bold green]OK All outputs written to:[/] {out_dir.resolve()}")
    rprint()

    # Print a summary of written files
    if out_dir.exists():
        for f in sorted(out_dir.iterdir()):
            if f.is_file():
                size_kb = f.stat().st_size / 1024
                rprint(f"    {f.name:<25s}  {size_kb:>8.1f} KB")

    rprint("\n[bold]Done![/]")
