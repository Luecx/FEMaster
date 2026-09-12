"""Legacy/non-element concentrated point-mass feature.

``PointMass`` applies mass, inertia and spring properties directly to a node
region through the native ``*POINTMASS`` feature.  It remains distinct from the
``MASS`` / ``ROTARYI`` / ``SPRING1`` point-element path because the solver input
semantics are different even when they can represent similar physical effects.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .feature import Feature


class PointMass(Feature):
    """Concentrated mass/inertia/spring feature on a named node region."""

    def __init__(
        self,
        node_region: str,
        mass: float = 0.0,
        *,
        inertia: Iterable[float] = (0.0, 0.0, 0.0),
        spring: Iterable[float] = (0.0, 0.0, 0.0),
        rotational_spring: Iterable[float] = (0.0, 0.0, 0.0),
    ) -> None:
        self.node_region = str(node_region)
        self.mass = float(mass)
        self.inertia = self._vec3(inertia, "inertia")
        self.spring = self._vec3(spring, "spring")
        self.rotational_spring = self._vec3(
            rotational_spring,
            "rotational_spring",
        )

    @staticmethod
    def _vec3(
        values: Iterable[float],
        name: str,
    ) -> tuple[float, float, float]:
        result = tuple(float(value) for value in values)
        if len(result) != 3:
            raise ValueError(f"{name} requires exactly 3 components")
        return result

    def export(self) -> str:
        return block([
            keyword("POINTMASS", NSET=self.node_region),
            csv((
                self.mass,
                *self.inertia,
                *self.spring,
                *self.rotational_spring,
            )),
        ])
