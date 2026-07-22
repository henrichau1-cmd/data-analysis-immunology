# Phase 0 + 1: Data Foundation and EDA — Implementation Plan

> **Note on format:** The writing-plans skill assumes an AI agent executes the
> plan. Here the *human author* is the implementer — that is the point of the
> project. Steps are therefore written as a curriculum. Steps marked
> **`YOU WRITE THIS`** are for the author to implement; unmarked code is
> boilerplate given outright so time is spent on ideas rather than typing.
> Claude's role is to review, unblock, and explain — not to write the marked code.

**Goal:** Build a tested, reproducible pipeline that fetches the MSK-IMPACT
immunotherapy cohort from the cBioPortal API, reshapes it into an analysis-ready
table, derives TMB features, and produces an exploratory notebook with figures.

**Architecture:** One-directional data flow — `fetch → clean → features → viz`.
All logic lives in importable modules under `src/`; notebooks contain narrative
and figures only. Raw API responses are cached to disk and never mutated.

**Tech Stack:** Python 3.13, pandas 2.3, numpy 2.3, matplotlib 3.10, seaborn
0.13, pytest 8.4, requests 2.32.

## Global Constraints

- Python 3.13.9 (Anaconda at `/opt/anaconda3`) — already installed
- pandas 2.3.3, numpy 2.3.5, matplotlib 3.10.6, seaborn 0.13.2, pytest 8.4.2,
  requests 2.32.5 — already installed; do not pin lower
- `lifelines` is NOT installed and is NOT needed in this plan (Phase 2 only)
- Study ID is exactly `tmb_mskcc_2018`
- API base is exactly `https://www.cbioportal.org/api` — no API key
- The S3 bulk download returns HTTP 403 and must not be used
- Expected cohort size: **1,661 patients, 1,661 samples, 11 cancer types**
- `data/raw/` and `data/interim/` are gitignored; `figures/` is committed
- Tests must run without network access
- Commit after every task; the git history is a project deliverable

---

## File Structure

| File | Responsibility |
| --- | --- |
| `src/data/fetch.py` | HTTP + caching only. Knows nothing about data meaning. |
| `src/data/clean.py` | Reshaping only. Long→wide pivot, dtype coercion, join. |
| `src/features/tmb.py` | TMB transformations as pure functions. |
| `src/viz/plots.py` | Figure styling and generation. |
| `tests/test_clean.py` | Pivot, parsing, and join tests. |
| `tests/test_tmb.py` | Feature transformation tests. |
| `tests/test_fetch.py` | Caching behavior tests (no network). |
| `notebooks/01-eda.ipynb` | Narrative and figures. No logic. |
| `Makefile` | `data`, `process`, `figures`, `test`, `all` |

---

# PHASE 0 — Acquisition and Environment

## Task 1: Repository scaffolding

**Files:**
- Create: `requirements.txt`, `Makefile`, `README.md`, `src/**/__init__.py`,
  `tests/__init__.py`, `data/.gitkeep`, `figures/.gitkeep`

- [ ] **Step 1: Create the directory tree**

Run these one line at a time. Each is self-contained — no backslash
continuations, which get mangled when pasted into a terminal and can silently
run the fragments as separate commands.

```bash
cd "/Users/henrichau/workspace/Data Analysis Immunology Project"
```

```bash
mkdir -p src/data src/features src/analysis src/models src/viz tests notebooks figures data/raw data/interim data/processed
```

```bash
touch src/__init__.py src/data/__init__.py src/features/__init__.py src/analysis/__init__.py src/models/__init__.py src/viz/__init__.py tests/__init__.py figures/.gitkeep data/processed/.gitkeep
```

Verify — this must print exactly **seven** paths (six packages under `src/`,
plus `tests/`):

```bash
find . -name "__init__.py" | sort
```

```
./src/__init__.py
./src/analysis/__init__.py
./src/data/__init__.py
./src/features/__init__.py
./src/models/__init__.py
./src/viz/__init__.py
./tests/__init__.py
```

- [ ] **Step 2: Write `requirements.txt`**

```
pandas>=2.3
numpy>=2.3
matplotlib>=3.10
seaborn>=0.13
scipy>=1.16
scikit-learn>=1.7
requests>=2.32
pytest>=8.4
jupyter
```

- [ ] **Step 3: Write `Makefile`**

Note: Makefile recipes require **tab** indentation, not spaces.

```makefile
.PHONY: data process figures test all clean

data:
	python -m src.data.fetch

process:
	python -m src.data.clean

figures:
	python -m src.viz.plots

test:
	pytest -v

all: data process figures

clean:
	rm -rf data/interim/* data/processed/*.csv figures/*.png
```

- [ ] **Step 4: Verify pytest runs on an empty suite**

Run: `pytest -v`
Expected: `no tests ran` — exit code 5. This confirms pytest is wired up.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "chore: scaffold project structure and build targets"
```

---

## Task 2: cBioPortal fetch with caching

**Files:**
- Create: `src/data/fetch.py`
- Test: `tests/test_fetch.py`

**Interfaces:**
- Produces: `fetch_clinical_data(data_type: str, cache_dir: Path = RAW_DIR,
  force: bool = False) -> list[dict]` where `data_type` is `"PATIENT"` or
  `"SAMPLE"`. Returns the parsed JSON list. Later tasks consume this output.

**Background — what the API returns.** Each element is a flat dict. Verified
live on 2026-07-22:

```python
# data_type="PATIENT"
{"uniquePatientKey": "UC0wMDAwMDU3...", "patientId": "P-0000057",
 "studyId": "tmb_mskcc_2018", "clinicalAttributeId": "OS_MONTHS", "value": "0"}

# data_type="SAMPLE"  — note it also carries patientId, which is the join key
{"uniqueSampleKey": "UC0wMDAwMDU3...", "uniquePatientKey": "UC0wMDAwMDU3...",
 "sampleId": "P-0000057-T01-IM3", "patientId": "P-0000057",
 "studyId": "tmb_mskcc_2018", "clinicalAttributeId": "GENE_PANEL",
 "value": "IMPACT341"}
```

One row per *attribute per patient*, not one row per patient. That is why Task 3
exists.

- [ ] **Step 1: Read the given test file**

Create `tests/test_fetch.py` exactly as below. These tests are given so you can
see what good tests look like before writing your own in Task 4.

```python
import json
from pathlib import Path

import pytest

from src.data import fetch


def test_uses_cache_and_makes_no_network_call(tmp_path, monkeypatch):
    """A cached file must be returned without touching the network."""
    cached = [{"patientId": "P-1", "clinicalAttributeId": "OS_MONTHS", "value": "5"}]
    (tmp_path / "PATIENT.json").write_text(json.dumps(cached))

    def explode(*args, **kwargs):
        raise AssertionError("network was called despite a valid cache")

    monkeypatch.setattr(fetch.requests, "get", explode)

    assert fetch.fetch_clinical_data("PATIENT", cache_dir=tmp_path) == cached


def test_force_bypasses_cache(tmp_path, monkeypatch):
    """force=True must re-download even when a cache exists."""
    (tmp_path / "PATIENT.json").write_text(json.dumps([{"stale": True}]))
    fresh = [{"fresh": True}]

    class FakeResponse:
        status_code = 200
        def raise_for_status(self): pass
        def json(self): return fresh

    monkeypatch.setattr(fetch.requests, "get", lambda *a, **k: FakeResponse())

    assert fetch.fetch_clinical_data("PATIENT", cache_dir=tmp_path, force=True) == fresh


def test_rejects_unknown_data_type(tmp_path):
    """Typos must fail loudly, not silently fetch nothing."""
    with pytest.raises(ValueError, match="PATIENT|SAMPLE"):
        fetch.fetch_clinical_data("PATIENTS", cache_dir=tmp_path)


def test_does_not_write_cache_on_failure(tmp_path, monkeypatch):
    """A failed request must leave no partial cache file behind.

    A truncated cache file is worse than no cache: the next run reads it,
    finds valid JSON, and silently analyses incomplete data.
    """
    class FakeResponse:
        status_code = 500
        def raise_for_status(self): raise RuntimeError("500 Server Error")
        def json(self): return None

    monkeypatch.setattr(fetch.requests, "get", lambda *a, **k: FakeResponse())
    monkeypatch.setattr(fetch.time, "sleep", lambda _: None)  # skip real backoff

    with pytest.raises(RuntimeError):
        fetch.fetch_clinical_data("PATIENT", cache_dir=tmp_path)

    assert not (tmp_path / "PATIENT.json").exists()


def test_retries_then_succeeds(tmp_path, monkeypatch):
    """One transient failure must not abort the run."""
    calls = {"n": 0}

    class FakeResponse:
        def __init__(self, ok): self.ok = ok
        def raise_for_status(self):
            if not self.ok: raise RuntimeError("503 Service Unavailable")
        def json(self): return [{"ok": True}]

    def flaky(*a, **k):
        calls["n"] += 1
        return FakeResponse(ok=calls["n"] > 1)      # fail once, then succeed

    monkeypatch.setattr(fetch.requests, "get", flaky)
    monkeypatch.setattr(fetch.time, "sleep", lambda _: None)

    assert fetch.fetch_clinical_data("PATIENT", cache_dir=tmp_path) == [{"ok": True}]
    assert calls["n"] == 2
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest tests/test_fetch.py -v`
Expected: all five FAIL with `ModuleNotFoundError` or `AttributeError` — the
module does not exist yet. **This step matters.** A test that has never failed
might be passing for the wrong reason.

- [ ] **Step 3: Create the module skeleton**

Create `src/data/fetch.py` with everything except the core function:

```python
"""Fetch clinical data from the cBioPortal REST API, with on-disk caching.

This module knows about HTTP and files. It knows nothing about what the data
means — no TMB, no survival, no cancer types. That lives downstream.
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
    "the endpoint is gone" - and only the second one should stop the pipeline.
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
    raise NotImplementedError("YOU WRITE THIS — see Step 4")


if __name__ == "__main__":
    for dt in VALID_TYPES:
        records = fetch_clinical_data(dt)
        print(f"{dt}: {len(records)} records cached")
```

- [ ] **Step 4: `YOU WRITE THIS` — implement `fetch_clinical_data`**

Replace the `raise NotImplementedError` line with your implementation. Roughly
10 lines. It must:

1. Reject a `data_type` not in `VALID_TYPES` with `ValueError`
2. Return the parsed contents of `cache_dir/<data_type>.json` if it exists and
   `force` is False
3. Otherwise call `_get_with_retry(_url(data_type))` — the given helper already
   handles retries, backoff, and raising on failure
4. Create `cache_dir` if needed, then write the JSON — **only after** a
   successful response
5. Return the records

**The decision that matters — step ordering.** Writing the cache file before
confirming the request succeeded produces a truncated file that *looks* valid on
the next run, and every downstream number is then silently wrong. The fourth
test exists to catch exactly this. Get the ordering right and the test passes for
the right reason.

**Hint:** `Path.exists()`, `Path.read_text()`, `Path.write_text()`,
`json.loads()`, `json.dumps()`, `Path.mkdir(parents=True, exist_ok=True)`.

- [ ] **Step 5: Run the tests**

Run: `pytest tests/test_fetch.py -v`
Expected: 5 passed.

- [ ] **Step 6: Fetch the real data**

Run: `make data`
Expected output:

```
PATIENT: 9966 records cached
SAMPLE: 24076 records cached
```

Both counts were verified against the live API on 2026-07-22.
The PATIENT count must be exactly **9966**. If it differs, stop and investigate
before continuing — the API may have changed and every downstream number depends
on this.

Verify the cache works:

```bash
ls -la data/raw/
make data     # must be instant the second time — no network
```

- [ ] **Step 7: Commit**

```bash
git add src/data/fetch.py tests/test_fetch.py
git commit -m "feat: add cBioPortal API client with on-disk caching"
```

---

# PHASE 1 — Wrangling and EDA

## Task 3: Long-to-wide pivot

**Files:**
- Create: `src/data/clean.py`
- Test: `tests/test_clean.py`

**Interfaces:**
- Consumes: `fetch.fetch_clinical_data` from Task 2
- Produces: `pivot_clinical(records: list[dict], id_field: str) -> pd.DataFrame`
  — one row per entity, one column per `clinicalAttributeId`, with `id_field`
  (`"patientId"` or `"sampleId"`) as a regular column, not the index.

**Background.** 9,966 records describe 1,661 patients — about 6 attributes each.
You need 1,661 rows and one column per attribute. This is the single most common
reshaping operation in real data work.

- [ ] **Step 1: Write the given tests**

Create `tests/test_clean.py`:

```python
import pandas as pd
import pytest

from src.data import clean

RECORDS = [
    {"patientId": "P-1", "clinicalAttributeId": "OS_MONTHS", "value": "5.5"},
    {"patientId": "P-1", "clinicalAttributeId": "SEX", "value": "Male"},
    {"patientId": "P-2", "clinicalAttributeId": "OS_MONTHS", "value": "12.0"},
    {"patientId": "P-2", "clinicalAttributeId": "SEX", "value": "Female"},
]


def test_pivot_produces_one_row_per_patient():
    out = clean.pivot_clinical(RECORDS, id_field="patientId")
    assert len(out) == 2


def test_pivot_creates_one_column_per_attribute():
    out = clean.pivot_clinical(RECORDS, id_field="patientId")
    assert {"OS_MONTHS", "SEX"}.issubset(out.columns)


def test_pivot_keeps_id_as_a_column_not_an_index():
    """Downstream merges need patientId as a real column."""
    out = clean.pivot_clinical(RECORDS, id_field="patientId")
    assert "patientId" in out.columns


def test_pivot_preserves_values():
    out = clean.pivot_clinical(RECORDS, id_field="patientId")
    row = out.loc[out["patientId"] == "P-1"].iloc[0]
    assert row["SEX"] == "Male"
    assert row["OS_MONTHS"] == "5.5"


def test_pivot_handles_missing_attribute_for_one_patient():
    """P-2 has no SEX. That must become NaN, not an error or a dropped row."""
    partial = [r for r in RECORDS if not (r["patientId"] == "P-2" and r["clinicalAttributeId"] == "SEX")]
    out = clean.pivot_clinical(partial, id_field="patientId")
    assert len(out) == 2
    assert pd.isna(out.loc[out["patientId"] == "P-2", "SEX"].iloc[0])
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_clean.py -v`
Expected: 5 FAIL, `ModuleNotFoundError: No module named 'src.data.clean'`.

- [ ] **Step 3: Create the module skeleton**

Create `src/data/clean.py`:

```python
"""Reshape raw cBioPortal records into analysis-ready tables.

This module changes the *shape* of data, not its meaning. It does not know what
TMB is or which direction survival runs.
"""
from pathlib import Path

import pandas as pd

from src.data.fetch import fetch_clinical_data

PROCESSED_DIR = Path(__file__).resolve().parents[2] / "data" / "processed"

EXPECTED_PATIENTS = 1661
EXPECTED_SAMPLES = 1661


def pivot_clinical(records, id_field):
    """Reshape long-format clinical records into one row per entity.

    Args:
        records: list of dicts with keys id_field, clinicalAttributeId, value.
        id_field: "patientId" or "sampleId".

    Returns:
        DataFrame with one row per entity and one column per attribute,
        with id_field as a regular column.
    """
    raise NotImplementedError("YOU WRITE THIS — see Step 4")
```

- [ ] **Step 4: `YOU WRITE THIS` — implement `pivot_clinical`**

About 4 lines. Build a DataFrame from `records`, then reshape it.

**Two approaches, and the choice is yours:**

- `df.pivot(index=..., columns=..., values=...)` — fails loudly on duplicate
  entries
- `df.pivot_table(..., aggfunc="first")` — silently keeps one of the duplicates

**The trade-off:** `pivot` raises if any `(patient, attribute)` pair appears
twice. `pivot_table` silently picks one. For a dataset you are seeing for the
first time, which behavior do you want? Consider that a duplicate here would mean
you misunderstand the data, and you would rather learn that now than discover it
in Phase 3.

Afterward you will need `.reset_index()` so `patientId` is a column, and
`.rename_axis(columns=None)` to drop the leftover axis name pandas leaves behind.

- [ ] **Step 5: Run the tests**

Run: `pytest tests/test_clean.py -v`
Expected: 5 passed.

- [ ] **Step 6: Commit**

```bash
git add src/data/clean.py tests/test_clean.py
git commit -m "feat: reshape long-format clinical records to wide"
```

---

## Task 4: Type coercion and survival status parsing

**Files:**
- Modify: `src/data/clean.py`
- Modify: `tests/test_clean.py`

**Interfaces:**
- Produces: `parse_os_status(value: str) -> bool` (True = deceased) and
  `coerce_types(df: pd.DataFrame) -> pd.DataFrame`

**Background.** Every value from the API is a string. `OS_MONTHS` is `"11"`, not
`11.0`. `OS_STATUS` is the composite string `"1:DECEASED"` or `"0:LIVING"`.
Verified distribution: 832 deceased, 829 living.

- [ ] **Step 1: `YOU WRITE THIS` — write the tests yourself**

This is your first time writing tests unaided. Add to `tests/test_clean.py` at
least five tests covering:

1. `parse_os_status("1:DECEASED")` returns `True`
2. `parse_os_status("0:LIVING")` returns `False`
3. An unrecognized value such as `"2:UNKNOWN"` — you decide the behavior, then
   test it
4. `coerce_types` converts `OS_MONTHS` to a numeric dtype
5. `coerce_types` converts `TMB_NONSYNONYMOUS` to a numeric dtype

**The decision that matters — what to do with an unrecognized status.**
Returning `False` treats an unknown patient as alive, which quietly inflates your
survival estimates and would never announce itself. Raising an exception stops
the pipeline on data you do not understand. Returning `None`/`NaN` propagates the
uncertainty honestly but forces every downstream consumer to handle it.

Pick one, write the test that pins it down, and add a one-line comment saying
why. There is no universally correct answer — but there is a correct *process*,
which is deciding deliberately rather than by accident.

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_clean.py -v`
Expected: your five new tests FAIL; the five from Task 3 still pass.

- [ ] **Step 3: `YOU WRITE THIS` — implement both functions**

Add to `src/data/clean.py`:

```python
NUMERIC_COLUMNS = ["OS_MONTHS", "TMB_NONSYNONYMOUS", "MUTATION_COUNT",
                   "AGE_AT_SEQ_REPORT"]


def parse_os_status(value):
    """Convert cBioPortal survival status to a boolean. True means deceased."""
    # YOU WRITE THIS


def coerce_types(df):
    """Return a copy of df with numeric columns converted from strings."""
    # YOU WRITE THIS
```

For `coerce_types`: iterate `NUMERIC_COLUMNS`, skip any not present, and use
`pd.to_numeric(df[col], errors="coerce")`. Work on a copy (`df = df.copy()`) so
you never mutate a caller's DataFrame — silent mutation is a genuinely nasty
class of bug.

- [ ] **Step 4: Run the tests**

Run: `pytest tests/test_clean.py -v`
Expected: 10 passed.

- [ ] **Step 5: Commit**

```bash
git add src/data/clean.py tests/test_clean.py
git commit -m "feat: add type coercion and survival status parsing"
```

---

## Task 5: Join patient and sample tables

**Files:**
- Modify: `src/data/clean.py`
- Modify: `tests/test_clean.py`

**Interfaces:**
- Produces: `build_analysis_table(force: bool = False) -> pd.DataFrame` — the
  merged, typed table written to `data/processed/analysis_table.csv`

**Background.** Patient-level attributes (`OS_MONTHS`, `OS_STATUS`, `DRUG_TYPE`,
`SEX`) and sample-level attributes (`TMB_NONSYNONYMOUS`, `CANCER_TYPE`,
`GENE_PANEL`) live in separate API responses. Both carry `patientId`, so the join
key is clean — but a join is still where row counts silently go wrong.

- [ ] **Step 1: Write the given regression tests**

Add to `tests/test_clean.py`:

```python
def test_merge_does_not_duplicate_rows():
    """One sample per patient in, one row out. Duplication here would silently
    weight some patients more heavily in every downstream statistic."""
    patients = pd.DataFrame({"patientId": ["P-1", "P-2"], "SEX": ["Male", "Female"]})
    samples = pd.DataFrame({"patientId": ["P-1", "P-2"],
                            "sampleId": ["P-1-T01", "P-2-T01"],
                            "CANCER_TYPE": ["Melanoma", "Glioma"]})
    out = clean.merge_patient_sample(patients, samples)
    assert len(out) == 2


def test_merge_does_not_silently_drop_patients():
    """An inner join on a mismatched key returns zero rows and looks like a
    legitimate empty result. Fail loudly instead."""
    patients = pd.DataFrame({"patientId": ["P-1"], "SEX": ["Male"]})
    samples = pd.DataFrame({"patientId": ["DIFFERENT-KEY"],
                            "sampleId": ["X-T01"], "CANCER_TYPE": ["Melanoma"]})
    with pytest.raises(ValueError):
        clean.merge_patient_sample(patients, samples)
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_clean.py -v`
Expected: the 2 new tests FAIL with `AttributeError: module has no attribute
'merge_patient_sample'`.

- [ ] **Step 3: `YOU WRITE THIS` — implement `merge_patient_sample`**

```python
def merge_patient_sample(patients, samples):
    """Merge sample-level onto patient-level data, one row per sample.

    Raises:
        ValueError: if the merge loses every row, or changes the row count.
    """
    # YOU WRITE THIS
```

About 8 lines: merge on `patientId`, then assert the result is non-empty and has
the same number of rows as `samples`. Raise `ValueError` with a message showing
the actual counts if not.

**Why the assertion belongs here rather than in a notebook:** this is precisely
the failure that destroyed the first dataset considered for this project — a
naive merge on mismatched keys returned zero rows and looked like a valid
result. An assertion inside the function protects every future caller, including
you in eight weeks when you have forgotten this detail.

- [ ] **Step 4: Run the tests**

Run: `pytest tests/test_clean.py -v`
Expected: 12 passed.

- [ ] **Step 5: Add the pipeline entry point**

Append to `src/data/clean.py`:

```python
def build_analysis_table(force=False):
    """Fetch, reshape, type, and merge into the analysis-ready table."""
    patients = pivot_clinical(fetch_clinical_data("PATIENT", force=force), "patientId")
    samples = pivot_clinical(fetch_clinical_data("SAMPLE", force=force), "sampleId")

    # Guard against silent schema drift: cBioPortal may change under us, and a
    # pipeline that quietly analyses 1,400 patients is worse than one that stops.
    if len(patients) != EXPECTED_PATIENTS:
        raise ValueError(f"expected {EXPECTED_PATIENTS} patients, got {len(patients)}")
    if len(samples) != EXPECTED_SAMPLES:
        raise ValueError(f"expected {EXPECTED_SAMPLES} samples, got {len(samples)}")

    merged = merge_patient_sample(patients, samples)
    merged = coerce_types(merged)
    merged["deceased"] = merged["OS_STATUS"].map(parse_os_status)

    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    merged.to_csv(PROCESSED_DIR / "analysis_table.csv", index=False)
    return merged


if __name__ == "__main__":
    table = build_analysis_table()
    print(f"rows={len(table)} cols={len(table.columns)}")
    print(f"cancer types={table['CANCER_TYPE'].nunique()}")
    print(f"deceased={int(table['deceased'].sum())}")
```

Note: the sample pivot uses `sampleId` as its id field, so `patientId` must
survive as a column for the merge. If it does not, revisit Task 3 — this is a
real wrinkle worth debugging rather than working around.

- [ ] **Step 6: Run the real pipeline**

Run: `make process`
Expected, matching the verified cohort:

```
rows=1661 cols=~24
cancer types=11
deceased=832
```

If `deceased` is not 832, your `parse_os_status` is wrong. Fix it before moving
on — every survival result in Phase 2 depends on this number.

- [ ] **Step 7: Commit**

```bash
git add src/data/clean.py tests/test_clean.py
git commit -m "feat: merge patient and sample tables into analysis table"
```

---

## Task 6: TMB features — the centerpiece

**Files:**
- Create: `src/features/tmb.py`
- Test: `tests/test_tmb.py`

**Interfaces:**
- Produces: `add_log_tmb(df) -> pd.DataFrame` adding column `log_tmb`;
  `add_within_cancer_percentile(df) -> pd.DataFrame` adding column `tmb_pct`;
  `add_tmb_high(df, threshold=0.8) -> pd.DataFrame` adding boolean `tmb_high`

**Background — why this task exists.** Verified TMB distribution: min 0.00,
median 5.87, p90 25.45, max 207.5. Severely right-skewed, so raw values are
unusable in most models and unreadable in most plots.

Worse, three sequencing panels of different sizes are in play — IMPACT341 (230
samples), IMPACT410 (1,001), IMPACT468 (430). A panel covering more genes finds
more mutations, so raw TMB is **not comparable across panels**. Samstein et al.
handle this by ranking TMB *within cancer type*, converting an
incomparable absolute value into a comparable relative one.

- [ ] **Step 1: Write the given tests**

Create `tests/test_tmb.py`:

```python
import numpy as np
import pandas as pd
import pytest

from src.features import tmb


def test_log_tmb_handles_zero():
    """TMB of exactly 0.00 is present in the real data. log(0) is -inf, which
    poisons every downstream mean and breaks most models."""
    df = pd.DataFrame({"TMB_NONSYNONYMOUS": [0.0, 5.0, 100.0]})
    out = tmb.add_log_tmb(df)
    assert np.isfinite(out["log_tmb"]).all()


def test_log_tmb_is_monotonic():
    """Order must be preserved: a bigger TMB must never get a smaller log_tmb."""
    df = pd.DataFrame({"TMB_NONSYNONYMOUS": [1.0, 5.0, 100.0]})
    out = tmb.add_log_tmb(df)
    assert out["log_tmb"].is_monotonic_increasing


def test_percentile_is_within_unit_interval():
    df = pd.DataFrame({
        "CANCER_TYPE": ["Melanoma"] * 5 + ["Glioma"] * 5,
        "TMB_NONSYNONYMOUS": [1.0, 2, 3, 4, 5, 10.0, 20, 30, 40, 50],
    })
    out = tmb.add_within_cancer_percentile(df)
    assert out["tmb_pct"].between(0, 1).all()


def test_percentile_is_computed_within_group_not_globally():
    """The whole point. Glioma values are all larger than Melanoma values, yet
    each group must span the full 0-1 range independently."""
    df = pd.DataFrame({
        "CANCER_TYPE": ["Melanoma"] * 5 + ["Glioma"] * 5,
        "TMB_NONSYNONYMOUS": [1.0, 2, 3, 4, 5, 10.0, 20, 30, 40, 50],
    })
    out = tmb.add_within_cancer_percentile(df)
    for cancer in ("Melanoma", "Glioma"):
        grp = out.loc[out["CANCER_TYPE"] == cancer, "tmb_pct"]
        assert grp.max() == pytest.approx(1.0)
        assert grp.min() == pytest.approx(0.2)


def test_percentile_does_not_change_row_count():
    df = pd.DataFrame({
        "CANCER_TYPE": ["Melanoma"] * 3,
        "TMB_NONSYNONYMOUS": [1.0, 2.0, 3.0],
    })
    assert len(tmb.add_within_cancer_percentile(df)) == 3


def test_tmb_high_flags_top_fraction():
    """threshold=0.8 means "top 20%". With 10 evenly ranked samples the
    percentiles are 0.1 ... 1.0, so exactly two samples (0.9, 1.0) qualify.
    The sample sitting exactly at 0.8 is the 8th of 10 - it is in the top 30%,
    not the top 20% - so the comparison must be strict."""
    df = pd.DataFrame({
        "CANCER_TYPE": ["Melanoma"] * 10,
        "TMB_NONSYNONYMOUS": [float(i) for i in range(1, 11)],
    })
    out = tmb.add_tmb_high(tmb.add_within_cancer_percentile(df), threshold=0.8)
    assert out["tmb_high"].sum() == 2
```

- [ ] **Step 2: Run to verify failure**

Run: `pytest tests/test_tmb.py -v`
Expected: 6 FAIL, `ModuleNotFoundError`.

- [ ] **Step 3: Create the skeleton**

Create `src/features/tmb.py`:

```python
"""Feature transformations for tumor mutational burden.

Every function takes a DataFrame and returns a NEW DataFrame with one column
added. No function mutates its input.
"""
import numpy as np
import pandas as pd


def add_log_tmb(df):
    """Add `log_tmb`, a log-transformed TMB that tolerates zeros."""
    # YOU WRITE THIS


def add_within_cancer_percentile(df):
    """Add `tmb_pct`: each sample's TMB percentile within its own cancer type."""
    # YOU WRITE THIS


def add_tmb_high(df, threshold=0.8):
    """Add boolean `tmb_high` for samples at or above the percentile threshold."""
    # YOU WRITE THIS
```

- [ ] **Step 4: `YOU WRITE THIS` — implement `add_log_tmb`**

Two lines. The decision: **how to handle TMB of exactly 0.**

- `np.log1p(x)` computes `log(1 + x)`, so zero maps to zero. Standard, and
  requires no arbitrary constant.
- `np.log(x + c)` for a small `c` you choose — but the choice of `c` is arbitrary
  and shifts your results.
- Dropping zeros — discards real patients and biases the cohort.

Pick one and leave a comment explaining the choice. Reviewers of a portfolio repo
notice this kind of comment.

- [ ] **Step 5: `YOU WRITE THIS` — implement `add_within_cancer_percentile`**

**This is the most important function in Phase 1.** Two lines, but they encode
the paper's core methodological insight.

You need, for each row, its TMB's percentile rank *among rows sharing its cancer
type*. The shape you want is:

```python
df.groupby("CANCER_TYPE")["TMB_NONSYNONYMOUS"].rank(pct=True)
```

`rank(pct=True)` returns rank as a fraction of group size. Because it is called
on a groupby, ranking happens inside each group and the result aligns back to the
original index — which is why row count is preserved and you can assign it
directly as a column.

**Why this beats a loop.** You could iterate the 11 cancer types, subset, rank,
and concatenate. That works, is ~10 lines, and is easy to get subtly wrong when
the index is not unique. The groupby version is one line and correct by
construction. Internalizing this pattern — *split-apply-combine* — is the main
skill of Phase 1.

Check `test_percentile_is_computed_within_group_not_globally` to confirm you
understand the expected result before implementing.

- [ ] **Step 6: `YOU WRITE THIS` — implement `add_tmb_high`**

One line. Note the test requires a **strict** `>` comparison, not `>=` — read the
docstring in `test_tmb_high_flags_top_fraction` for why. Boundary conditions like
this are exactly what tests exist to pin down; the difference is invisible in a
plot but changes which patients enter your Phase 2 comparison.

- [ ] **Step 7: Run the tests**

Run: `pytest tests/test_tmb.py -v`
Expected: 6 passed.

- [ ] **Step 8: Sanity-check against real data**

```bash
python -c "
import pandas as pd
from src.features import tmb
df = pd.read_csv('data/processed/analysis_table.csv')
df = tmb.add_tmb_high(tmb.add_within_cancer_percentile(tmb.add_log_tmb(df)))
print(df.groupby('CANCER_TYPE')['tmb_high'].agg(['sum','size','mean']).round(3))
"
```

Expected: every cancer type shows a `mean` near 0.2, because you flagged the top
20% *within each type*. If one type shows 0.9 and another 0.02, you ranked
globally rather than within groups.

- [ ] **Step 9: Commit**

```bash
git add src/features/tmb.py tests/test_tmb.py
git commit -m "feat: add TMB log transform and within-cancer-type percentile"
```

---

## Task 7: Vectorization exercise

**Files:**
- Create: `notebooks/00-vectorization.ipynb`

This task produces no library code. It exists to make one lesson concrete.

- [ ] **Step 1: Write the loop version**

In a notebook cell, compute z-scored TMB within each cancer type using an
explicit Python loop over `df.iterrows()`. Time it with `%%timeit`.

- [ ] **Step 2: `YOU WRITE THIS` — write the vectorized version**

Same computation using `groupby().transform()`. Roughly:

```python
g = df.groupby("CANCER_TYPE")["TMB_NONSYNONYMOUS"]
df["tmb_z"] = (df["TMB_NONSYNONYMOUS"] - g.transform("mean")) / g.transform("std")
```

Time it the same way.

- [ ] **Step 3: Verify both produce identical results**

```python
assert np.allclose(loop_result, vectorized_result, equal_nan=True)
```

Getting the same answer two ways is how you build confidence in the fast version.

- [ ] **Step 4: Write a markdown cell recording the speedup**

Expect somewhere between 100x and 1000x on 1,661 rows. Note the ratio and one
sentence on why: the loop runs Python bytecode per row, while the vectorized
version dispatches to compiled C once per group.

- [ ] **Step 5: Commit**

```bash
git add notebooks/00-vectorization.ipynb
git commit -m "docs: add vectorization benchmark notebook"
```

---

## Task 8: Plotting module and EDA notebook

**Files:**
- Create: `src/viz/plots.py`, `notebooks/01-eda.ipynb`
- Output: `figures/*.png`

**Interfaces:**
- Consumes: `data/processed/analysis_table.csv`, all three functions from Task 6
- Produces: `set_style() -> None` and one function per figure, each returning a
  matplotlib `Figure`

- [ ] **Step 1: Create the styling module**

Create `src/viz/plots.py`:

```python
"""Figure generation. Consumes results; produces PNGs. No analysis logic."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")          # render without a display, so `make figures` works
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from src.features import tmb

FIGURES_DIR = Path(__file__).resolve().parents[2] / "figures"
PROCESSED = Path(__file__).resolve().parents[2] / "data" / "processed" / "analysis_table.csv"


def set_style():
    sns.set_theme(style="whitegrid", context="talk")
    plt.rcParams["figure.dpi"] = 150
    plt.rcParams["savefig.bbox"] = "tight"


def load_features():
    df = pd.read_csv(PROCESSED)
    return tmb.add_tmb_high(tmb.add_within_cancer_percentile(tmb.add_log_tmb(df)))
```

- [ ] **Step 2: `YOU WRITE THIS` — the skew figure**

Write `plot_tmb_distribution(df) -> plt.Figure` producing two side-by-side
panels: raw TMB on the left, `log_tmb` on the right.

This figure is the argument for the log transform. Someone should look at it and
immediately understand why raw TMB is unusable — the left panel should be a spike
at zero with a long invisible tail out to 207.5.

- [ ] **Step 3: `YOU WRITE THIS` — TMB by cancer type**

Write `plot_tmb_by_cancer_type(df) -> plt.Figure`. A boxplot or violin of
`log_tmb` across the 11 cancer types.

**A real design decision:** sort the cancer types by median TMB rather than
alphabetically. Alphabetical ordering carries no information; ordering by the
value being compared lets a reader rank groups at a glance. Small choices like
this separate a figure that communicates from one that merely displays.

Expect melanoma and NSCLC high (UV and smoking mutagenesis), glioma low. If your
figure shows that, your pipeline is working — it reproduces known biology.

- [ ] **Step 4: `YOU WRITE THIS` — cohort composition**

Write `plot_cohort_overview(df) -> plt.Figure` showing sample counts by cancer
type and by drug class. Verified totals to check against: NSCLC 350, Melanoma
320, Bladder 215; PD-1/PD-L1 1,307, Combo 255, CTLA-4 99.

- [ ] **Step 5: Add the figure-writing entry point**

```python
def main():
    set_style()
    df = load_features()
    FIGURES_DIR.mkdir(exist_ok=True)
    for name, fn in [
        ("tmb_distribution", plot_tmb_distribution),
        ("tmb_by_cancer_type", plot_tmb_by_cancer_type),
        ("cohort_overview", plot_cohort_overview),
    ]:
        fig = fn(df)
        fig.savefig(FIGURES_DIR / f"{name}.png")
        plt.close(fig)
        print(f"wrote figures/{name}.png")


if __name__ == "__main__":
    main()
```

- [ ] **Step 6: Generate figures**

Run: `make figures`
Expected: three PNGs written to `figures/`. Open them and check they are readable
— axis labels present, no overlapping tick text, log axis clearly labelled as log.

- [ ] **Step 7: `YOU WRITE THIS` — build the EDA notebook**

Create `notebooks/01-eda.ipynb`. It imports from `src/` and contains **no
transformation logic**. Structure:

1. Question being asked
2. Cohort description — how many patients, which cancers, which drugs
3. Each figure, with a markdown cell interpreting it
4. At least 3 more figures of your own choosing — survival status by cancer type,
   TMB by drug class, age distribution, panel-versus-TMB, whatever you find
   interesting
5. A closing section: what you learned, and what you now want to test statistically
   in Phase 2

Target 6–8 figures total. Every figure needs a sentence saying what it shows;
a figure without interpretation is decoration.

- [ ] **Step 8: Commit**

```bash
git add src/viz/plots.py notebooks/01-eda.ipynb figures/*.png
git commit -m "feat: add plotting module and exploratory analysis notebook"
```

---

## Task 9: README and phase closeout

**Files:**
- Modify: `README.md`

- [ ] **Step 1: `YOU WRITE THIS` — write the README**

Sections, per the spec:

1. Question and headline observation, with `figures/tmb_by_cancer_type.png`
   embedded near the top
2. Data provenance — Samstein et al., *Nature Genetics* 2019, accessed via the
   cBioPortal REST API, study `tmb_mskcc_2018`, 1,661 patients
3. What exists so far (Phase 1 complete: pipeline and EDA; Phases 2–3 pending)
4. Reproduction: clone, `pip install -r requirements.txt`, `make all`
5. Repository layout
6. Limitations — retrospective cohort, no control arm, panel heterogeneity,
   no direct immune measurements

Write section 6 honestly. It is the section most likely to impress someone
evaluating you.

- [ ] **Step 2: Verify reproducibility from scratch**

```bash
mv data/raw data/raw_backup
make all
```

Everything must regenerate from nothing. This is the claim your README makes; verify
it before making it.

```bash
rm -rf data/raw_backup
```

- [ ] **Step 3: Run the full suite**

Run: `pytest -v`
Expected: 23 passed (5 fetch + 12 clean + 6 tmb), no network access required.

- [ ] **Step 4: Commit and push**

```bash
git add -A
git commit -m "docs: add README with results and reproduction instructions"
git remote add origin <your-github-url>
git push -u origin main
```

---

## Definition of Done

- [ ] `pytest -v` → 23 passed, no network needed
- [ ] `make all` regenerates every figure from an empty `data/`
- [ ] `data/processed/analysis_table.csv` has 1,661 rows, 11 cancer types, 832 deceased
- [ ] Within-cancer-type `tmb_high` rate is ~0.2 in every cancer type
- [ ] `notebooks/01-eda.ipynb` has 6–8 figures, each interpreted
- [ ] README present with an honest limitations section
- [ ] Repository pushed publicly with incremental commit history
- [ ] You can explain every line you wrote without referring to this plan

## What You Should Be Able to Do Afterward

- Reshape long ↔ wide without looking it up
- Reach for `groupby().transform()` when a per-group value must align back to rows
- Explain why vectorized operations beat loops, with a measured number
- Write a test that pins down a decision rather than restating the implementation
- Recognize that a merge returning zero rows is a bug, not a result
