"""Format-independent FEMaster solution, frame and field hierarchy.

Result files are read directly through ``Result.read_res`` / ``Result.read_frd``;
there is intentionally no separate public I/O package.
"""

from .result import Result
from .result_frame import Frame
from .result_solution import Solution

__all__ = ["Frame", "Result", "Solution"]
