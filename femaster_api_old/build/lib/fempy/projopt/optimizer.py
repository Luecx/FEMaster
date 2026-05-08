from __future__ import annotations
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import List, Optional, Tuple, Dict, Any
import numpy as np
from .types import Array, Real
from .objective import Objective
from .constraint import Constraint, Relation

try:
    from scipy.optimize import minimize, Bounds
except Exception as e:
    minimize = None
    Bounds = None

class Method(Enum):
    TRUST_CONSTR = auto()
    SLSQP        = auto()
    L_BFGS_B     = auto()  # bounds only

    def to_scipy(self) -> str:
        if self is Method.TRUST_CONSTR: return "trust-constr"
        if self is Method.SLSQP:        return "SLSQP"
        if self is Method.L_BFGS_B:     return "L-BFGS-B"
        raise ValueError("Unknown method")

@dataclass
class OptParams:
    # Choose solver
    method: Method = Method.TRUST_CONSTR
    # Common knobs (we map the relevant ones to SciPy depending on method)
    maxiter: int = 500
    gtol: float = 1e-8     # for trust-constr
    ftol: float = 1e-12    # for SLSQP / L-BFGS-B
    verbose_solver: int = 0  # passed to trust-constr 'verbose'
    disp: bool = False       # passed to SLSQP/L-BFGS-B 'disp'
    # Bounds:
    bounds: Tuple[Array, Array] | None = None    # explicit arrays
    lb_val: float | None = None                  # uniform lower bound (applied at solve-time using x0 size)
    ub_val: float | None = None                  # uniform upper bound
    # Printing a compact iteration table from our wrapper
    verbose: bool = True

@dataclass
class OptStats:
    success: bool
    status: int
    message: str
    niter: int
    nfev: int | None
    njev: int | None
    nhev: int | None
    final_grad_inf: float
    final_fun: float

@dataclass
class IterInfo:
    k: int
    f: Real
    step: Real
    grad_inf: Real

class SciPyOptimizer:
    """
    Minimal optimizer: wraps scipy.optimize.minimize, using your Objective/Constraint layers.
    Captures per-iteration histories for plotting and exposes stats via get_stats().
    """
    def __init__(self, params: OptParams = OptParams()):
        if minimize is None:
            raise RuntimeError("scipy is required. Please `pip install scipy`.")
        self.P = params
        self.objective: Objective | None = None
        self.constraints: List[Constraint] = []
        self.responses = []  # optional logging in the print table
        # histories (per SciPy iteration)
        self.history_f: List[float] = []
        self.history_gradinf: List[float] = []
        self.history_step: List[float] = []
        self._prev_x: Array | None = None
        self._last_result: Any | None = None  # scipy OptimizeResult

    def set_objective(self, obj: Objective) -> None:
        self.objective = obj

    def add_constraint(self, c: Constraint) -> None:
        self.constraints.append(c)

    def add_response(self, r) -> None:
        self.responses.append(r)

    def _build_bounds(self, n: int):
        # Priority: explicit arrays > (lb_val/ub_val) > None
        if self.P.bounds is not None:
            lb, ub = self.P.bounds
            return Bounds(np.asarray(lb, float), np.asarray(ub, float))
        if self.P.lb_val is not None or self.P.ub_val is not None:
            lb = np.full(n, -np.inf, dtype=float) if self.P.lb_val is None else np.full(n, float(self.P.lb_val))
            ub = np.full(n, +np.inf, dtype=float) if self.P.ub_val is None else np.full(n, float(self.P.ub_val))
            return Bounds(lb, ub)
        return None

    def solve(self, x: Array) -> IterInfo:
        assert self.objective and self.objective.valid(), "Objective not set."
        x0 = np.asarray(x, dtype=float)

        def fun(xv):
            return float(self.objective.value(xv))

        def jac(xv):
            g = self.objective.grad(xv)
            return np.asarray(g, dtype=float)

        # SciPy constraint format:
        #  - 'ineq' means fun(x) >= 0 is FEASIBLE
        # Our internal inequalities are LT (g(x) <= 0 feasible) => pass -g(x) as 'ineq'
        scipy_constraints = []
        for c in self.constraints:
            if c.rel == Relation.EQ:
                scipy_constraints.append({
                    "type": "eq",
                    "fun":  lambda xv, cf=c: float(cf.g(xv)),
                    "jac":  lambda xv, cf=c: np.asarray(cf.grad(xv), dtype=float)
                })
            else:
                scipy_constraints.append({
                    "type": "ineq",
                    "fun":  lambda xv, cf=c: float(-cf.g(xv)),
                    "jac":  lambda xv, cf=c: np.asarray(-cf.grad(xv), dtype=float)
                })

        bounds = self._build_bounds(x0.size)

        if self.P.verbose:
            self._print_header()
            f0 = fun(x0)
            g0 = jac(x0)
            self._record(0, x0, f0, float(np.linalg.norm(g0, ord=np.inf)), step=np.nan)

        k_counter = {"k": 0}
        def cb(xk):
            k_counter["k"] += 1
            f = fun(xk)
            g = jac(xk)
            grad_inf = float(np.linalg.norm(g, ord=np.inf))
            step = np.nan
            if self._prev_x is not None:
                step = float(np.linalg.norm(xk - self._prev_x))
            self._prev_x = xk.copy()
            self._record(k_counter["k"], xk, f, grad_inf, step)

        method = self.P.method.to_scipy()
        options: Dict[str, Any] = {}
        if method == "trust-constr":
            options.update({"maxiter": self.P.maxiter, "gtol": self.P.gtol, "verbose": self.P.verbose_solver})
        elif method in ("SLSQP", "L-BFGS-B"):
            options.update({"maxiter": self.P.maxiter, "ftol": self.P.ftol, "disp": self.P.disp})
        else:
            options.update({"maxiter": self.P.maxiter})

        res = minimize(fun, x0, method=method, jac=jac,
                       constraints=scipy_constraints,
                       bounds=bounds,
                       options=options,
                       callback=cb if method in ("trust-constr", "SLSQP") else None)

        self._last_result = res

        # Final record (ensure last is captured)
        k_final = k_counter["k"]
        f_final = float(res.fun)
        g_final = jac(res.x)
        grad_inf = float(np.linalg.norm(g_final, ord=np.inf))
        if len(self.history_f) == 0 or self.history_f[-1] != f_final:
            step = np.nan
            if self._prev_x is not None:
                step = float(np.linalg.norm(res.x - self._prev_x))
            self._record(k_final + 1, res.x, f_final, grad_inf, step)

        if self.P.verbose:
            self._print_row(k_final + 1, f_final, np.nan, grad_inf, res.x)

        x[:] = res.x  # write back
        nit = int(getattr(res, "nit", k_final))
        return IterInfo(k=nit, f=f_final, step=np.nan, grad_inf=grad_inf)

    def get_stats(self) -> OptStats:
        """Return final stats collected from SciPy's OptimizeResult."""
        if self._last_result is None:
            raise RuntimeError("solve() has not been called yet.")
        res = self._last_result
        nfev = getattr(res, "nfev", None)
        njev = getattr(res, "njev", None)
        nhev = getattr(res, "nhev", None)
        g = None
        try:
            g = self.objective.grad(res.x)
            grad_inf = float(np.linalg.norm(g, ord=np.inf))
        except Exception:
            grad_inf = float("nan")
        return OptStats(
            success=bool(getattr(res, "success", False)),
            status=int(getattr(res, "status", 0)),
            message=str(getattr(res, "message", "")),
            niter=int(getattr(res, "nit", 0)),
            nfev=int(nfev) if nfev is not None else None,
            njev=int(njev) if njev is not None else None,
            nhev=int(nhev) if nhev is not None else None,
            final_grad_inf=grad_inf,
            final_fun=float(getattr(res, "fun", np.nan)),
        )

    # ---- printing / history ----
    def _record(self, k: int, x: Array, f: float, grad_inf: float, step: float):
        self.history_f.append(float(f))
        self.history_gradinf.append(float(grad_inf))
        self.history_step.append(float(step))
        if self.P.verbose:
            self._print_row(k, f, step, grad_inf, x)

    def _print_header(self):
        head = f"{'it':>4}  {'f':>14}  {'step':>10}  {'||grad||_inf':>14}"
        for r in self.responses:
            head += f"  {r.name:>14}"
        print(head)

    def _print_row(self, k: int, f: Real, step: Real, grad_inf: Real, x: Array):
        row = f"{k:4d}  {f:14.6e}  {step:10.3e}  {grad_inf:14.6e}"
        for r in self.responses:
            row += f"  {r.value(x):14.6e}"
        print(row)
