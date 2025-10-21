from __future__ import annotations
import numpy as np
from typing import Callable

Array = np.ndarray
Real  = float
FunF  = Callable[[Array], Real]   # scalar function f(x)
GradF = Callable[[Array], Array]  # gradient ∇f(x)
