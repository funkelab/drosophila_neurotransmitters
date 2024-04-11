from enum import IntEnum, Enum
import pandas as pd
from pydantic import BaseModel


class Ternary(IntEnum):
    false = -1
    unknown = 0
    true = 1


class SourceType(Enum):
    prediction = "prediction"
    unknown = "unknown"
    easi_fish = "easi_fish"
    immunohistochemistry = "immunohistochemistry"
    other = "other"


class Evidence(IntEnum):
    """
    How strong the evidence is for the presence of a neurotransmitter in a cell type. 
    """
    limited = 1
    moderate = 2"
    strong = 3


class Row(BaseModel):
    cell_type: str
    acetylcholine: Ternary
    glutamate: Ternary
    gaba: Ternary
    dopamine: Ternary
    serotonin: Ternary
    octopamine: Ternary
    tyramine: Ternary
    glycine: Ternary
    histamine: Ternary
    nitric_oxide: Ternary
    source: str
    source_type: SourceType
    evidence: Evidence


def validate_row(row: dict) -> Row:
    return Row(**row)


def validate_df(df):
    return df.apply(validate_row, axis=1)


def test_example():
    filename = "data.csv"
    df = pd.read_csv(filename)
    validate_df(df)
