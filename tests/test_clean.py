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

def test_parse_os_status_dtype_changed():
    df = pd.DataFrame({"OS_MONTHS": ["5.5", "12.0"]})
    out = clean.coerce_types(df)
    assert pd.api.types.is_numeric_dtype(out["OS_MONTHS"])

def test_parse_os_status_deceased_is_true():
    assert clean.parse_os_status("1:DECEASED") is True

def test_parse_os_status_not_deceased_is_false():
    assert clean.parse_os_status("0:LIVING") is False

def test_parse_os_status_unknown_is_ValueError():
    with pytest.raises(ValueError):
        clean.parse_os_status("2:UNKNOWN")

def test_coerce_types_makes_tmb_numeric():
    df = pd.DataFrame({"TMB_NONSYNONYMOUS": ["5.87", "25.45"]})
    out = clean.coerce_types(df)
    assert pd.api.types.is_numeric_dtype(out["TMB_NONSYNONYMOUS"])
