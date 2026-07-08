"""
tcga/api.py — cBioPortal REST API v2 client
============================================

Thin wrappers around the public cBioPortal REST API (no authentication
required for public studies).

Provides:
  - list_studies()                   List all available studies
  - get_study(study_id)              Validate a study ID and return metadata
  - get_molecular_profiles(id)       List molecular profiles for a study
  - detect_available_data_types(id)  Check which data types exist for a study
  - get_clinical_data(id, type)      Fetch sample- or patient-level clinical data
  - get_sample_list(id)              Get sample IDs for a study
  - get_molecular_data(profile_id, sample_ids, genes)
                                     Fetch molecular data for specific profile
  - get_mutations(study_id, profile_id, sample_ids)
                                     Fetch mutation data
"""

from __future__ import annotations

import time
from typing import Any, Optional

import requests

BASE_URL = "https://www.cbioportal.org/api"
_SESSION = requests.Session()
_SESSION.headers.update({"Accept": "application/json"})


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get(
    endpoint: str,
    params: Optional[dict[str, Any]] = None,
    retries: int = 3,
    timeout: int = 60,
) -> Any:
    """
    Make a GET request to the cBioPortal API with simple retry logic.
    Raises requests.HTTPError on failure after all retries.
    """
    url = f"{BASE_URL}/{endpoint.lstrip('/')}"
    last_exc: Optional[Exception] = None

    for attempt in range(retries):
        try:
            response = _SESSION.get(url, params=params, timeout=timeout)
            response.raise_for_status()
            return response.json()
        except (requests.ConnectionError, requests.Timeout) as exc:
            last_exc = exc
            if attempt < retries - 1:
                time.sleep(2 ** attempt)  # exponential back-off
        except requests.HTTPError:
            raise

    raise last_exc  # type: ignore[misc]


def _post(
    endpoint: str,
    json_body: Any = None,
    params: Optional[dict[str, Any]] = None,
    retries: int = 3,
    timeout: int = 120,
) -> Any:
    """
    Make a POST request to the cBioPortal API with simple retry logic.
    The API uses POST for fetching large payloads (molecular data).
    """
    url = f"{BASE_URL}/{endpoint.lstrip('/')}"
    last_exc: Optional[Exception] = None

    for attempt in range(retries):
        try:
            response = _SESSION.post(url, json=json_body, params=params, timeout=timeout)
            response.raise_for_status()
            return response.json()
        except (requests.ConnectionError, requests.Timeout) as exc:
            last_exc = exc
            if attempt < retries - 1:
                time.sleep(2 ** attempt)
        except requests.HTTPError:
            raise

    raise last_exc  # type: ignore[misc]


# ---------------------------------------------------------------------------
# Public API functions
# ---------------------------------------------------------------------------

def list_studies() -> list[dict]:
    """
    Return all publicly available studies as a list of dicts.

    Each dict has keys:
      studyId, name, description, cancerType, allSampleCount
    """
    raw = _get(
        "/studies",
        params={"projection": "SUMMARY", "pageSize": 10_000, "pageNumber": 0},
    )
    return [
        {
            "studyId": s["studyId"],
            "name": s.get("name", ""),
            "description": s.get("description", ""),
            "cancerType": s.get("cancerType", {}).get("name", "")
                if isinstance(s.get("cancerType"), dict)
                else "",
            "allSampleCount": s.get("allSampleCount", 0),
        }
        for s in raw
    ]


def get_study(study_id: str) -> Optional[dict]:
    """
    Return metadata for a specific study, or None if the study ID is not found.
    """
    try:
        return _get(f"/studies/{study_id}")
    except requests.HTTPError as exc:
        if exc.response is not None and exc.response.status_code == 404:
            return None
        raise


def get_molecular_profiles(study_id: str) -> list[dict]:
    """
    Return a list of molecular profiles (data types) available for a study.

    Each profile has keys such as molecularProfileId, name, molecularAlterationType.
    """
    return _get(f"/studies/{study_id}/molecular-profiles")


# ---------------------------------------------------------------------------
# Data type detection
# ---------------------------------------------------------------------------

# Maps our canonical data type names to cBioPortal molecularAlterationType values
_ALTERATION_TYPE_MAP: dict[str, list[str]] = {
    "rnaseq":    ["MRNA_EXPRESSION"],
    "cna":       ["COPY_NUMBER_ALTERATION"],
    "mutations": ["MUTATION_EXTENDED"],
    "rppa":      ["PROTEIN_LEVEL"],
}

# Preferred profile selection: when multiple profiles exist for a type,
# pick the one whose ID contains one of these keywords (in priority order).
_PROFILE_PREFERENCES: dict[str, list[str]] = {
    "rnaseq": ["rna_seq_v2_mrna"],      # raw RSEM, not z-scores
    "cna":    ["gistic"],                 # discrete (-2 to +2), not log2
    "rppa":   ["rppa_Zscores", "rppa"],   # z-scores preferred
}


def detect_available_data_types(study_id: str) -> dict[str, bool]:
    """
    Inspect a study's molecular profiles to determine which data types are
    available.  Clinical data is assumed to always be present.

    Returns a dict: {data_type: bool}
    """
    try:
        profiles = get_molecular_profiles(study_id)
    except requests.HTTPError:
        return {key: True for key in ["clinical", *_ALTERATION_TYPE_MAP]}

    profile_types = {p.get("molecularAlterationType", "") for p in profiles}

    availability: dict[str, bool] = {"clinical": True}
    for data_type, alt_types in _ALTERATION_TYPE_MAP.items():
        availability[data_type] = any(at in profile_types for at in alt_types)

    return availability


def resolve_profile_id(study_id: str, data_type: str) -> Optional[str]:
    """
    Find the best molecular profile ID for a given canonical data type.

    For example, for 'rnaseq' on study 'kirc_tcga_pan_can_atlas_2018',
    this returns 'kirc_tcga_pan_can_atlas_2018_rna_seq_v2_mrna'.
    """
    profiles = get_molecular_profiles(study_id)
    alt_types = _ALTERATION_TYPE_MAP.get(data_type, [])

    # Filter to matching alteration types
    candidates = [
        p for p in profiles
        if p.get("molecularAlterationType", "") in alt_types
    ]

    if not candidates:
        return None

    # If there are preferences, pick the first matching preference
    prefs = _PROFILE_PREFERENCES.get(data_type, [])
    for pref in prefs:
        for c in candidates:
            pid = c["molecularProfileId"]
            # For rnaseq, make sure we get the raw counts, not z-scores
            if data_type == "rnaseq" and "zscores" in pid.lower():
                continue
            if pref.lower() in pid.lower():
                return pid

    # Fall back: for rnaseq, filter out z-score profiles
    if data_type == "rnaseq":
        non_zscore = [c for c in candidates if "zscore" not in c["molecularProfileId"].lower()]
        if non_zscore:
            return non_zscore[0]["molecularProfileId"]

    # Final fall back: first candidate
    return candidates[0]["molecularProfileId"]


# ---------------------------------------------------------------------------
# Data fetching
# ---------------------------------------------------------------------------

def get_sample_ids(study_id: str) -> list[str]:
    """Return all sample IDs for a study."""
    samples = _get(f"/studies/{study_id}/samples", params={"projection": "ID", "pageSize": 10_000})
    return [s["sampleId"] for s in samples]


def get_clinical_data(
    study_id: str,
    clinical_data_type: str = "SAMPLE",
) -> list[dict]:
    """
    Fetch clinical data (SAMPLE or PATIENT level) for a study.

    Returns a list of dicts with keys:
      sampleId/patientId, clinicalAttributeId, value
    """
    records = _get(
        f"/studies/{study_id}/clinical-data",
        params={
            "clinicalDataType": clinical_data_type,
            "projection": "DETAILED",
            "pageSize": 100_000,
            "pageNumber": 0,
        },
    )
    return records


def get_molecular_data(
    molecular_profile_id: str,
    sample_ids: list[str],
    entry_gene_ids: Optional[list[int]] = None,
) -> list[dict]:
    """
    Fetch molecular data (expression, CNA, RPPA) for a given profile.

    Uses the POST /molecular-profiles/{id}/molecular-data/fetch endpoint
    to handle large sample lists. Requests are chunked to avoid 502 Bad Gateway
    timeouts on the public API.
    """
    from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn

    all_records = []
    chunk_size = 100  # Smaller chunks to avoid timeouts
    total_samples = len(sample_ids)

    with Progress(
        TextColumn("    [blue]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("Fetching chunks", total=total_samples)
        for i in range(0, total_samples, chunk_size):
            chunk = sample_ids[i:i + chunk_size]
            body: dict[str, Any] = {"sampleIds": chunk}
            if entry_gene_ids:
                body["entrezGeneIds"] = entry_gene_ids

            records = _post(
                f"/molecular-profiles/{molecular_profile_id}/molecular-data/fetch",
                json_body=body,
                params={"projection": "SUMMARY"},
            )
            all_records.extend(records)
            progress.advance(task, len(chunk))

    return all_records


def get_mutations(
    molecular_profile_id: str,
    sample_ids: list[str],
) -> list[dict]:
    """
    Fetch mutation data for a given mutation profile.

    Uses the POST /molecular-profiles/{id}/mutations/fetch endpoint.
    Returns one record per mutation call. Requests are chunked.
    """
    from rich.progress import Progress, BarColumn, TextColumn, TimeElapsedColumn

    all_records = []
    chunk_size = 100
    total_samples = len(sample_ids)

    with Progress(
        TextColumn("    [blue]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("Fetching chunks", total=total_samples)
        for i in range(0, total_samples, chunk_size):
            chunk = sample_ids[i:i + chunk_size]
            body = {"sampleIds": chunk}
            records = _post(
                f"/molecular-profiles/{molecular_profile_id}/mutations/fetch",
                json_body=body,
                params={"projection": "DETAILED", "pageSize": 100_000, "pageNumber": 0},
            )
            all_records.extend(records)
            progress.advance(task, len(chunk))
        
    return all_records


def get_all_genes() -> list[dict]:
    """Return a list of all genes known to cBioPortal (for resolving entrezGeneId → Hugo_Symbol)."""
    return _get("/genes", params={"projection": "SUMMARY", "pageSize": 100_000})
