"""Backend execution and result-reading API."""

from .femaster import FEMaster, ProcessResult
from .result_reader import Frame, LoadCase, Result, ResultReader

__all__ = [
    "FEMaster",
    "Frame",
    "LoadCase",
    "ProcessResult",
    "Result",
    "ResultReader",
]
