"""Orthotropic engineering-constant elasticity."""

from ..common.format import block, csv, keyword
from .elasticity import Elasticity


class OrthotropicElasticity(Elasticity):
    """Orthotropic engineering constants in material coordinates."""

    def __init__(
        self,
        e1: float,
        e2: float,
        e3: float,
        nu12: float,
        nu13: float,
        nu23: float,
        g12: float,
        g13: float,
        g23: float,
    ) -> None:
        self.values = tuple(float(value) for value in (
            e1, e2, e3, nu12, nu13, nu23, g12, g13, g23
        ))

    def export(self) -> str:
        return block([
            keyword("ELASTIC", TYPE="ENGINEERINGCONSTANTS"),
            csv(self.values),
        ])
