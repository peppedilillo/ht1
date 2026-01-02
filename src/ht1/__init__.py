"""
fm1trig - Transient search for SpIRIT/HERMES-FM1 SRA data.
"""

from ht1.ht1 import (
    # Enums
    Quadrant,
    EnBand,
    ErrorCode,
    TriggerStatus,
    # Type aliases
    BandData,
    Data,
    Hit,
    Index,
    Interval,
    # Specs and exceptions
    SRASpecs,
    InvalidSRA,
    # Functions
    sra_parse,
    ma_range,
    moving_average,
    significance,
    search_qbdata,
    search_data,
    search_filepath,
    hit_tointerval,
    summarize,
    # Classes
    TriggerDyadic,
)

__all__ = [
    "Quadrant",
    "EnBand",
    "ErrorCode",
    "TriggerStatus",
    "BandData",
    "Data",
    "Hit",
    "Index",
    "Interval",
    "SRASpecs",
    "InvalidSRA",
    "sra_parse",
    "ma_range",
    "moving_average",
    "significance",
    "search_qbdata",
    "search_data",
    "search_filepath",
    "hit_tointerval",
    "summarize",
    "TriggerDyadic",
]
