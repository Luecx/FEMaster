"""Legacy POINTMASS feature."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .feature import Feature


class PointMass(Feature):
    """Concentrated mass, inertia and spring properties on a node region."""

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
        self.inertia = self._vec3(inertia)
        self.spring = self._vec3(spring)
        self.rotational_spring = self._vec3(rotational_spring)

    def export(self) -> str:
        return block([
            keyword("POINTMASS", NSET=self.node_region),
            csv((self.mass, *self.inertia, *self.spring, *self.rotational_spring)),
        ])

    @staticmethod
    def _vec3(values: Iterable[float]) -> tuple[float, float, float]:
        result = tuple(float(value) for value in values)
        if len(result) != 3:
            raise ValueError("expected exactly 3 components")
        return result
