from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, Sequence
import numpy as np
from .types import Array, Real, FunF, GradF

@dataclass(frozen=True)
class Response:
    """
    A named scalar function r(x) with its gradient ∇r(x).
    Closed under linear and common nonlinear operations.
    """
    name: str
    value: FunF         # r(x)  -> scalar
    grad:  GradF        # ∇r(x) -> vector

    # ---------- linear ops ----------
    def __neg__(self) -> "Response":
        return Response(
            name=f"(-{self.name})",
            value=lambda x, v=self.value: -v(x),
            grad =lambda x, g=self.grad:  -g(x),
        )

    def __add__(self, other: "Response") -> "Response":
        return Response(
            name=f"({self.name}+{other.name})",
            value=lambda x, v1=self.value, v2=other.value: v1(x) + v2(x),
            grad =lambda x, g1=self.grad,  g2=other.grad : g1(x) + g2(x),
        )
    __radd__ = __add__

    def __sub__(self, other: "Response") -> "Response":
        return Response(
            name=f"({self.name}-{other.name})",
            value=lambda x, v1=self.value, v2=other.value: v1(x) - v2(x),
            grad =lambda x, g1=self.grad,  g2=other.grad : g1(x) - g2(x),
        )
    def __rsub__(self, other: "Response") -> "Response":
        return (-self).__add__(other)

    def __mul__(self, c: Real | "Response") -> "Response":
        # scalar * response OR response * response (product rule)
        if isinstance(c, (int, float)):
            c = float(c)
            return Response(
                name=f"({c}*{self.name})",
                value=lambda x, v=self.value, c=c: c * v(x),
                grad =lambda x, g=self.grad,  c=c: c * g(x),
            )
        elif isinstance(c, Response):
            return Response(
                name=f"({self.name}*{c.name})",
                value=lambda x, v1=self.value, v2=c.value: v1(x) * v2(x),
                grad =lambda x, g1=self.grad,  g2=c.grad, v1=self.value, v2=c.value:
                v2(x) * g1(x) + v1(x) * g2(x),
            )
        else:
            raise TypeError("Unsupported type for multiplication")
    __rmul__ = __mul__

    def __truediv__(self, other: Real | "Response") -> "Response":
        eps = 1e-16
        if isinstance(other, (int, float)):
            c = float(other)
            if abs(c) < eps: raise ZeroDivisionError("division by ~0 scalar")
            inv = 1.0 / c
            return Response(
                name=f"({self.name}/{c})",
                value=lambda x, v=self.value, inv=inv: inv * v(x),
                grad =lambda x, g=self.grad,  inv=inv: inv * g(x),
            )
        elif isinstance(other, Response):
            # quotient rule: (u/v)' = (v*u' - u*v')/v^2
            def val(x, u=self.value, v=other.value):
                vx = v(x)
                if abs(vx) < eps:
                    vx = eps if vx >= 0 else -eps
                return u(x) / vx
            def grd(x, gu=self.grad, gv=other.grad, u=self.value, v=other.value):
                ux = u(x); vx = v(x)
                if abs(vx) < eps:
                    vx = eps if vx >= 0 else -eps
                return (vx * gu(x) - ux * gv(x)) / (vx * vx)
            return Response(name=f"({self.name}/{other.name})", value=val, grad=grd)
        else:
            raise TypeError("Unsupported type for division")

    # power by scalar: (r**p)' = p*r^{p−1}*r'
    def __pow__(self, p: Real) -> "Response":
        p = float(p)
        def val(x, v=self.value, p=p):
            rx = v(x)
            return rx ** p
        def grd(x, g=self.grad, v=self.value, p=p):
            rx = v(x)
            if p == 0.0:
                return 0.0 * g(x)
            return p * (rx ** (p - 1.0)) * g(x)
        return Response(name=f"({self.name}**{p})", value=val, grad=grd)

    # ---------- nonlinear composition helpers ----------
    def map(self, name: str, phi: Callable[[Real], Real], dphi: Callable[[Real], Real]) -> "Response":
        """Unary nonlinear map φ(r(x)) with chain rule."""
        return Response(
            name=f"{name}({self.name})",
            value=lambda x, v=self.value, phi=phi: phi(v(x)),
            grad =lambda x, g=self.grad,  v=self.value, dphi=dphi: dphi(v(x)) * g(x),
        )

    @staticmethod
    def combine(name: str,
                responses: Sequence["Response"],
                phi: Callable[[Sequence[Real]], Real],
                dphi: Callable[[Sequence[Real]], Sequence[Real]]) -> "Response":
        """Multi-arg composition y = φ(r1,...,rk) with provided partials ∂φ/∂r_i."""
        rs = list(responses)
        def val(x):
            vals = [r.value(x) for r in rs]
            return float(phi(vals))
        def grd(x):
            vals = [r.value(x) for r in rs]
            partials = dphi(vals)
            gsum = np.zeros_like(rs[0].grad(x))
            for wi, ri in zip(partials, rs):
                gsum += wi * ri.grad(x)
            return gsum
        return Response(name=name, value=val, grad=grd)

def compose_response(name: str, value: FunF, grad: GradF) -> Response:
    return Response(name=name, value=value, grad=grad)
