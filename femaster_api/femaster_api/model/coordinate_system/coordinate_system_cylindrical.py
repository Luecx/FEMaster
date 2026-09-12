"""Cylindrical coordinate system with explicit origin, axis and reference.

All three defining vectors are required and validated as three-dimensional.
Keeping the cylindrical definition in its own class prevents rectangular and
cylindrical parameter conventions from being mixed in one generic orientation
object.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .coordinate_system import CoordinateSystem


class CylindricalCoordinateSystem(CoordinateSystem):
    """Named cylindrical FEMaster coordinate system."""

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

    @staticmethod
    def _vec3(values: Iterable[float]) -> tuple[float, float, float]:
        result = tuple(float(value) for value in values)
        if len(result) != 3:
            raise ValueError("expected exactly 3 vector components")
        return result

    def export(self) -> str:
        return block([
            keyword("ORIENTATION", NAME=self.name, TYPE="CYLINDRICAL"),
            csv((*self.origin, *self.axis, *self.reference)),
        ])
