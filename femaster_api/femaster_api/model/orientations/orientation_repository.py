"""Repository for coordinate-system definitions."""

from __future__ import annotations

from femaster_api.utils import normalize_name

from .cylindrical_orientation import CylindricalOrientation
from .rectangular_orientation import RectangularOrientation

Vector3 = tuple[float, float, float]
Orientation = RectangularOrientation | CylindricalOrientation


class OrientationRepository:
    """Repository for named coordinate systems."""

    def __init__(self) -> None:
        self._items: dict[str, Orientation] = {}

    def add(self, orientation: Orientation) -> Orientation:
        if not isinstance(orientation, (RectangularOrientation, CylindricalOrientation)):
            raise TypeError("orientation must be a RectangularOrientation or CylindricalOrientation")
        self._items[normalize_name(orientation.name)] = orientation
        return orientation

    def get(self, name: str) -> Orientation:
        key = normalize_name(name)
        try:
            return self._items[key]
        except KeyError as exc:
            raise KeyError(f"unknown orientation: {key}") from exc

    def __getitem__(self, name: str) -> Orientation:
        return self.get(name)

    def has(self, value: str | Orientation) -> bool:
        if isinstance(value, (RectangularOrientation, CylindricalOrientation)):
            return any(item is value for item in self._items.values())
        return normalize_name(value) in self._items

    def all(self) -> tuple[Orientation, ...]:
        return tuple(self._items[key] for key in sorted(self._items))

    def __iter__(self):
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._items)
