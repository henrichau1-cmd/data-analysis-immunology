"""Fetch clinical data from the cBioPortal REST API, with on-disk caching.

This module knows about HTTP and files. It knows nothing about what the data
means -- no TMB, no survival, no cancer types. That lives downstream.
"""
import json
import time
from pathlib import Path

import requests

API_BASE = "https://www.cbioportal.org/api"
STUDY_ID = "tmb_mskcc_2018"
RAW_DIR = Path(__file__).resolve().parents[2] / "data" / "raw"
VALID_TYPES = ("PATIENT", "SAMPLE")
PAGE_SIZE = 200000
TIMEOUT = 120
MAX_ATTEMPTS = 3


def _url(data_type: str) -> str:
    return (
        f"{API_BASE}/studies/{STUDY_ID}/clinical-data"
        f"?clinicalDataType={data_type}&pageSize={PAGE_SIZE}"
    )


def _get_with_retry(url):
    """GET with exponential backoff. Raises RuntimeError if all attempts fail.

    Transient network failures are common on a 24k-record request. Retrying
    three times with a growing pause distinguishes "the wifi blipped" from
    "the endpoint is gone" -- and only the second one should stop the pipeline.
    """
    last_error = None
    for attempt in range(MAX_ATTEMPTS):
        try:
            response = requests.get(url, timeout=TIMEOUT)
            response.raise_for_status()
            return response.json()
        except Exception as error:          # noqa: BLE001 - deliberately broad
            last_error = error
            if attempt < MAX_ATTEMPTS - 1:
                time.sleep(2 ** attempt)     # 1s, then 2s
    raise RuntimeError(f"failed after {MAX_ATTEMPTS} attempts: {last_error}")


def fetch_clinical_data(data_type, cache_dir=RAW_DIR, force=False):
    """Return clinical data for the study, using a cached copy when available.

    Args:
        data_type: "PATIENT" or "SAMPLE".
        cache_dir: directory holding <data_type>.json.
        force: if True, ignore any cache and re-download.

    Returns:
        list[dict] of clinical records.

    Raises:
        ValueError: if data_type is not PATIENT or SAMPLE.
        RuntimeError: if the HTTP request fails.
    """
    if data_type not in VALID_TYPES:
        raise ValueError("data_type must be PATIENT or SAMPLE")
    path = Path(cache_dir) / f"{data_type}.json"
    if path.exists() and not force:
        return json.loads(path.read_text())

    records = _get_with_retry(_url(data_type))
    Path(cache_dir).mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(records))
    return records

if __name__ == "__main__":
    for dt in VALID_TYPES:
        records = fetch_clinical_data(dt)
        print(f"{dt}: {len(records)} records cached")
