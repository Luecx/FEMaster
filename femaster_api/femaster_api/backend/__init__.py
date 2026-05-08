"""FEMaster backend execution and result-reading API."""

from .result_reader import FEMasterResult, Frame, LoadCaseResult, LoadcaseResult, Result, ResultField, ResultReader
from .runner import FEMaster, FEMasterRunner, RunResult, find_femaster_executable

__all__ = [
    "FEMaster",
    "FEMasterResult",
    "FEMasterRunner",
    "Frame",
    "LoadCaseResult",
    "LoadcaseResult",
    "Result",
    "ResultField",
    "ResultReader",
    "RunResult",
    "find_femaster_executable",
]
