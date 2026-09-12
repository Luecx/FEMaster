"""Public entry point of the FEMaster Python API.

The package re-exports the semantic model hierarchy from ``femaster_api.model``.
There is intentionally no separate public ``io`` layer: editable input models
are read through ``Project.read_inp`` and result files through
``Result.read_res`` / ``Result.read_frd``.
"""

from .model import *

__all__ = [name for name in globals() if not name.startswith("_")]
