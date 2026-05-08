from __future__ import annotations
from dataclasses import dataclass
from enum import Enum, auto
from .types import FunF, GradF, Real
from .responses import Response

class Relation(Enum):
    LT = auto()  # user: r(x) <= rhs  -> g(x)=r-rhs <=0
    GT = auto()  # user: r(x) >= rhs  -> g(x)=rhs-r <=0 (normalized)
    EQ = auto()  # user: r(x) == rhs  -> g(x)=r-rhs ==0

@dataclass
class Constraint:
    g:    FunF
    grad: GradF
    rel:  Relation
    name: str = ""

# Kept for API compatibility (not used by SciPy backend).
@dataclass
class ActiveSetRule:
    tol_eq:   Real = 1e-10
    tol_ineq: Real = 1e-8
    directional_activation: bool = True

def constraint_from_response(resp: Response, rel: Relation, rhs: Real, name: str | None = None) -> Constraint:
    rhs = float(rhs)
    if rel == Relation.LT:
        gfun   = lambda x: resp.value(x) - rhs
        gradfn = lambda x: resp.grad(x)
        out_rel = Relation.LT
    elif rel == Relation.GT:
        gfun   = lambda x: rhs - resp.value(x)
        gradfn = lambda x: -resp.grad(x)
        out_rel = Relation.LT   # normalize to LT form internally
    elif rel == Relation.EQ:
        gfun   = lambda x: resp.value(x) - rhs
        gradfn = lambda x: resp.grad(x)
        out_rel = Relation.EQ
    else:
        raise ValueError("unknown relation")

    cname = name if name else f"{resp.name}_{rel.name}_{rhs}"
    return Constraint(g=gfun, grad=gradfn, rel=out_rel, name=cname)
