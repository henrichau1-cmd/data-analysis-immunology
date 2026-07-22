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
