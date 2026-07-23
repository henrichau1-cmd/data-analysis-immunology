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
    pivoted = df.pivot(index=id_field, columns= "clinicalAttributeId", values= "value").reset_index()
    return pivoted
    # ------------------------------------------------------------------
    # YOU WRITE THIS -- about 4 lines. See Task 3, Step 4 in the plan.
    #
    # 1. Build a DataFrame from the list of dicts
    # 2. Reshape it: one row per id_field, one column per
    #    clinicalAttributeId, filled with value
    # 3. The reshape leaves id_field as the index -- make it a column again
    # 4. Return the result
    #
    # Useful: pd.DataFrame(), DataFrame.pivot(index=, columns=, values=),
    #         DataFrame.reset_index()
    # ------------------------------------------------------------------
    raise NotImplementedError("YOU WRITE THIS -- see Task 3, Step 4")
