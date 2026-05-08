"""Top-level façade for the FEMaster Python API."""

from ._version import __version__
from .backend.result_reader import Result, ResultReader
from .backend.runner import FEMaster, FEMasterRunner
from .export.femaster_writer import FEMasterWriter
from .model import Model

__all__ = [
    "__version__",
    "FEMaster",
    "FEMasterRunner",
    "FEMasterWriter",
    "Model",
    "Result",
    "ResultReader",
]
