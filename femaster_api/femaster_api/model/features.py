"""Assembly-level non-topological FEMaster features."""

from __future__ import annotations

from typing import Iterable, Iterator

from .._format import block, csv, keyword


def _vec3(values: Iterable[float]) -> tuple[float, float, float]:
    result = tuple(float(value) for value in values)
    if len(result) != 3:
        raise ValueError("expected exactly 3 components")
    return result


class Feature:
    """Base class for one assembly-level FEMaster feature."""

    def export(self) -> str:
        raise NotImplementedError


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
        self.node_region       = str(node_region)
        self.mass              = float(mass)
        self.inertia           = _vec3(inertia)
        self.spring            = _vec3(spring)
        self.rotational_spring = _vec3(rotational_spring)

    def export(self) -> str:
        return block([
            keyword("POINTMASS", NSET=self.node_region),
            csv((self.mass, *self.inertia, *self.spring, *self.rotational_spring)),
        ])


class FeatureRepository:
    """Ordered collection of heterogeneous non-topological features."""

    def __init__(self) -> None:
        self._items: list[Feature] = []

    def add(self, feature: Feature) -> Feature:
        if not isinstance(feature, Feature):
            raise TypeError("feature must derive from Feature")
        self._items.append(feature)
        return feature

    def export(self) -> str:
        return "\n\n".join(feature.export() for feature in self._items)

    def __iter__(self) -> Iterator[Feature]:
        return iter(self._items)

    def __len__(self) -> int:
        return len(self._items)
