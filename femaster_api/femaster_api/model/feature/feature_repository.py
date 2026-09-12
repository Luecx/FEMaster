"""Ordered heterogeneous repository of non-topological features.

Features currently do not require persistent names, so repository order is the
only ownership metadata.  Each feature owns its native serialization.
"""

from __future__ import annotations

from collections.abc import Iterator

from .feature import Feature


class FeatureRepository:
    """Own non-topological model features in insertion order."""

    def __init__(self) -> None:
        self._items: list[Feature] = []

    def add(self, feature: Feature) -> Feature:
        """Append a concrete feature and return it."""

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
