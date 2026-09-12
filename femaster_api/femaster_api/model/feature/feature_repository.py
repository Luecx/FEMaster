"""Ordered heterogeneous feature repository."""

from collections.abc import Iterator

from .feature import Feature


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
