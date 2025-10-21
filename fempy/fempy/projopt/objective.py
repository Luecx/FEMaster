from __future__ import annotations
from dataclasses import dataclass
from enum import Enum, auto
from typing import Sequence
import numpy as np
from .types import FunF, GradF, Real
from .responses import Response

@dataclass
class Objective:
    value: FunF
    grad:  GradF
    def valid(self) -> bool:
        return callable(self.value) and callable(self.grad)

class Sense(Enum):
    MIN = auto()
    MAX = auto()

def objective_from_response(resp: Response, sense: Sense = Sense.MIN) -> Objective:
    if sense == Sense.MIN:
        return Objective(value=resp.value, grad=resp.grad)
    else:
        return Objective(value=lambda x: -resp.value(x),
                         grad =lambda x: -resp.grad(x))

def objective_from_weighted_responses(resps: Sequence[Response],
                                      weights: Sequence[Real],
                                      sense: Sense = Sense.MIN) -> Objective:
    w = np.asarray(weights, dtype=float)
    assert len(resps) == len(w) and len(w) > 0
    def f(x): return float(sum(w[i] * resps[i].value(x) for i in range(len(resps))))
    def g(x):
        gsum = np.zeros_like(resps[0].grad(x))
        for i, r in enumerate(resps):
            gsum += w[i]*r.grad(x)
        return gsum
    if sense == Sense.MIN:
        return Objective(value=f, grad=g)
    else:
        return Objective(value=lambda x: -f(x), grad=lambda x: -g(x))
