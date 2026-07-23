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
    df = pd.DataFrame(records)
    pivoted = df.pivot(index=id_field, columns="clinicalAttributeId", values="value").reset_index()
    return pivoted


NUMERIC_COLUMNS = ["OS_MONTHS", "TMB_NONSYNONYMOUS", "MUTATION_COUNT",
                   "AGE_AT_SEQ_REPORT"]


def parse_os_status(value):
    """Convert cBioPortal survival status to a boolean. True means deceased."""
    if value == "1:DECEASED":
        return True
    elif value == "0:LIVING":
        return False
    else:
        raise ValueError("Unrecognized OS: {}".format(value))

def coerce_types(df):
    """Return a copy of df with numeric columns converted from strings."""
    df = df.copy()
    for column in NUMERIC_COLUMNS:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")
    return df
    # YOU WRITE THIS -- see Task 4, Step 3
    raise NotImplementedError("YOU WRITE THIS -- see Task 4, Step 3")
