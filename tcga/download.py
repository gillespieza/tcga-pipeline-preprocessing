"""
tcga/download.py — Study tarball download and extraction
=========================================================

Downloads a cBioPortal study data package (the same zip available on the
Datasets page), extracts it, and locates the relevant flat data files.

Caching:  if the extracted directory already exists, download is skipped
          unless force=True is passed.

Public API:
  download_study(study_id, cache_dir, force)  -> Path to extracted directory
  find_data_files(extracted_dir)              -> dict mapping type -> Path|None
"""

from __future__ import annotations

import re
import shutil
import tarfile
import zipfile
from pathlib import Path
from typing import Optional

import requests
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)

# cBioPortal study data download endpoint (same as "Data" button on website)
_DOWNLOAD_URL = (
    "https://www.cbioportal.org/study/downloadStudyData.do?studyId={study_id}"
)

# ---------------------------------------------------------------------------
# File pattern registry
# Each key maps to a compiled regex that matches the relevant cBioPortal file.
# Patterns are intentionally broad to handle naming variants across studies.
# ---------------------------------------------------------------------------
_FILE_PATTERNS: dict[str, re.Pattern] = {
    # Clinical
    "clinical_patient": re.compile(r"data_clinical_patient", re.IGNORECASE),
    "clinical_sample":  re.compile(r"data_clinical_sample",  re.IGNORECASE),
    # Genomic / molecular
    "rnaseq":    re.compile(r"data_mrna_seq",                  re.IGNORECASE),
    "cna":       re.compile(r"data_cna(?!.*seg)(?!.*log2)",    re.IGNORECASE),
    "mutations": re.compile(r"data_mutations",                  re.IGNORECASE),
    "rppa":      re.compile(r"data_rppa",                       re.IGNORECASE),
}


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------

def download_study(
    study_id: str,
    cache_dir: Path,
    force: bool = False,
) -> Path:
    """
    Download and extract a cBioPortal study data package.

    Parameters
    ----------
    study_id  : cBioPortal study identifier (e.g. 'kirc_tcga_pan_can_atlas_2018')
    cache_dir : Directory where downloads and extracted files are cached.
    force     : If True, delete any existing cached extraction and re-download.

    Returns
    -------
    Path to the extracted study directory containing raw data files.
    """
    extracted_dir = cache_dir / study_id
    archive_path  = cache_dir / f"{study_id}.zip"

    # Return cached extraction if available (and force not requested)
    if extracted_dir.exists() and any(extracted_dir.iterdir()) and not force:
        return extracted_dir

    # Clean up stale cache if forcing
    if force and extracted_dir.exists():
        shutil.rmtree(extracted_dir)

    extracted_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    url = _DOWNLOAD_URL.format(study_id=study_id)
    _stream_download(url, archive_path)
    _extract_archive(archive_path, extracted_dir)

    # Clean up the archive — keep only the extracted files
    archive_path.unlink(missing_ok=True)

    return extracted_dir


def find_data_files(extracted_dir: Path) -> dict[str, Optional[Path]]:
    """
    Walk the extracted study directory and match files to data type keys using
    filename pattern matching.

    Returns a dict: {data_type_key: Path or None}
    When multiple files match (rare), the shallowest one is preferred.
    """
    # Gather all .txt and .tsv files in the extracted tree
    candidates: list[Path] = sorted(
        list(extracted_dir.rglob("*.txt")) + list(extracted_dir.rglob("*.tsv")),
        key=lambda p: (len(p.parts), p.name),  # shallower first
    )

    result: dict[str, Optional[Path]] = {}
    for data_type, pattern in _FILE_PATTERNS.items():
        matches = [f for f in candidates if pattern.search(f.name)]
        result[data_type] = matches[0] if matches else None

    return result


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _stream_download(url: str, dest: Path) -> None:
    """Stream a file from url to dest with a rich progress bar."""
    with requests.get(url, stream=True, timeout=120, allow_redirects=True) as resp:
        resp.raise_for_status()
        total_bytes = int(resp.headers.get("content-length", 0)) or None

        with Progress(
            TextColumn("  [bold blue]{task.description}"),
            BarColumn(),
            DownloadColumn(),
            TransferSpeedColumn(),
            TimeRemainingColumn(),
        ) as progress:
            task = progress.add_task("Downloading study archive", total=total_bytes)

            with open(dest, "wb") as fh:
                for chunk in resp.iter_content(chunk_size=65_536):
                    fh.write(chunk)
                    progress.advance(task, len(chunk))


def _extract_archive(archive_path: Path, dest_dir: Path) -> None:
    """
    Extract a zip or tar.gz/tgz archive into dest_dir.
    Raises ValueError for unrecognised formats.
    """
    if zipfile.is_zipfile(archive_path):
        with zipfile.ZipFile(archive_path, "r") as zf:
            zf.extractall(dest_dir)
    else:
        try:
            with tarfile.open(archive_path) as tf:
                tf.extractall(dest_dir)
        except tarfile.TarError as exc:
            raise ValueError(
                f"Cannot extract '{archive_path.name}': unrecognised archive format. "
                f"Original error: {exc}"
            ) from exc
