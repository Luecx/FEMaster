"""Legacy/non-element concentrated point-mass feature on a ``NodeRegion``.

``PointMass`` applies mass, inertia and spring properties directly to an actual
node-region object.  Its native ``NSET`` token is derived from ``node_region.name``
only while exporting ``*POINTMASS``.  The feature therefore cannot silently
retain a dangling node-set name.

It remains distinct from element-based mass/inertia/spring properties because
the solver input semantics are different.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..region.region_node import NodeRegion
from .feature import Feature


class PointMass(Feature):
    """Concentrated mass/inertia/spring feature on one node region."""

    def __init__(
        self,
        node_region: NodeRegion,
        mass: float = 0.0,
        *,
        inertia: Iterable[float] = (0.0, 0.0, 0.0),
        spring: Iterable[float] = (0.0, 0.0, 0.0),
        rotational_spring: Iterable[float] = (0.0, 0.0, 0.0),
    ) -> None:
        if not isinstance(node_region, NodeRegion):
            raise TypeError("node_region must be a NodeRegion object")
        self.node_region = node_region
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
            keyword("POINTMASS", NSET=self.node_region.name),
            csv((
                self.mass,
                *self.inertia,
                *self.spring,
                *self.rotational_spring,
            )),
        ])
