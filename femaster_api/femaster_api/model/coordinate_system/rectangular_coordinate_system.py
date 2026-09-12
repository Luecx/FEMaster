"""Rectangular coordinate system."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .coordinate_system import CoordinateSystem


class RectangularCoordinateSystem(CoordinateSystem):
    """Rectangular system defined by one to three basis directions."""

    def __init__(
        self,
        name: str,
        x_axis: Iterable[float],
        y_axis: Iterable[float] | None = None,
        z_axis: Iterable[float] | None = None,
    ) -> None:
        super().__init__(name)
        self.x_axis = self._vec3(x_axis)
        self.y_axis = None if y_axis is None else self._vec3(y_axis)
        self.z_axis = None if z_axis is None else self._vec3(z_axis)

    def export(self) -> str:
        values: list[float] = list(self.x_axis)
        if self.y_axis is not None:
            values.extend(self.y_axis)
        if self.z_axis is not None:
            values.extend(self.z_axis)
        return block([
            keyword("ORIENTATION", NAME=self.name, TYPE="RECTANGULAR"),
            csv(values),
        ])

    @staticmethod
    def _vec3(values: Iterable[float]) -> tuple[float, float, float]:
        result = tuple(float(value) for value in values)
        if len(result) != 3:
            raise ValueError("expected exactly 3 vector components")
        return result
