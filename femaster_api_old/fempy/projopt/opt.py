from __future__ import annotations
import numpy as np
from typing import List
from .types import Array, Real, FunF, GradF
from .responses import Response, compose_response
from .objective import Objective, Sense, objective_from_response, objective_from_weighted_responses
from .constraint import Constraint, Relation, constraint_from_response
from .optimizer import SciPyOptimizer, OptParams, Method

# Module-level registry to auto-attach created responses (for logging table)
_RESP_REG: List[Response] = []

class Opt:
    # --- response creation (auto-register) ---
    @staticmethod
    def create_response(name: str, value: FunF, grad: GradF) -> Response:
        r = compose_response(name, value, grad)
        _RESP_REG.append(r)
        return r

    # --- objective builders ---
    @staticmethod
    def objective(resp_or_obj, sense: Sense = Sense.MIN) -> Objective:
        if isinstance(resp_or_obj, Response):
            return objective_from_response(resp_or_obj, sense)
        elif isinstance(resp_or_obj, Objective):
            return resp_or_obj
        else:
            raise TypeError("objective() expects a Response or Objective")

    @staticmethod
    def objective_weighted(responses, weights, sense: Sense = Sense.MIN) -> Objective:
        return objective_from_weighted_responses(responses, weights, sense)

    # --- constraint builder ---
    @staticmethod
    def constraint(resp_or_g, relation: Relation, rhs: Real | None = None, name: str | None = None) -> Constraint:
        if isinstance(resp_or_g, Response):
            assert rhs is not None, "rhs required for response-based constraint"
            return constraint_from_response(resp_or_g, relation, rhs, name)
        elif isinstance(resp_or_g, tuple) and len(resp_or_g) == 2:
            gfun, grad = resp_or_g
            if relation == Relation.GT:
                G = lambda x: -gfun(x)
                dG = lambda x: -grad(x)
                return Constraint(g=G, grad=dG, rel=Relation.LT, name=name or "nonlinear_GT")
            else:
                return Constraint(g=gfun, grad=grad, rel=relation, name=name or "nonlinear")
        else:
            raise TypeError("constraint() expects Response or (g, grad) tuple")

    # --- optimizer factory ---
    @staticmethod
    def new(params: OptParams = OptParams()) -> SciPyOptimizer:
        opt = SciPyOptimizer(params)
        for r in _RESP_REG:
            opt.add_response(r)
        return opt

    @staticmethod
    def clear_registered_responses():
        _RESP_REG.clear()
