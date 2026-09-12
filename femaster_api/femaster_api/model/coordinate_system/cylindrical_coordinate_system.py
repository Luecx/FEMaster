"""Cylindrical coordinate system."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .coordinate_system import CoordinateSystem


class CylindricalCoordinateSystem(CoordinateSystem):
    """Cylindrical system defined by origin, axis and reference direction."""

    def __init__(
        self,
        name: str,
        origin: Iterable[float],
        axis: Iterable[float],
        reference: Iterable[float],
    ) -> None:
        super().__init__(name)
        self.origin = self._vec3(origin)
        self.axis = self._vec3(axis)
        self.reference = self._vec3(reference)

    def export(self) -> str:
        return block([
            keyword("ORIENTATION", NAME=self.name, TYPE="CYLINDRICAL"),
            csv((*self.origin, *self.axis, *self.reference)),
        ])

    @staticmethod
    def _vec3(values: Iterable[float]) -> tuple[float, float, float]:
        result = tuple(float(value) for value in values)
        if len(result) != 3:
            raise ValueError("expected exactly 3 vector components")
        return result
