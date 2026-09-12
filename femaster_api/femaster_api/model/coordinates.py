"""Named global coordinate-system definitions."""

from __future__ import annotations

from typing import Iterable

from .._format import block, csv, keyword
from ..repository import NamedObject, NamedRepository


def _vec3(values: Iterable[float]) -> tuple[float, float, float]:
    result = tuple(float(value) for value in values)
    if len(result) != 3:
        raise ValueError("expected exactly 3 vector components")
    return result


class CoordinateSystem(NamedObject):
    """Base class for named FEMaster coordinate systems."""

    def to_femaster(self) -> str:
        raise NotImplementedError


class RectangularCoordinateSystem(CoordinateSystem):
    """Rectangular coordinate system defined by one to three basis directions."""

    def __init__(
        self,
        name: str,
        x_axis: Iterable[float],
        y_axis: Iterable[float] | None = None,
        z_axis: Iterable[float] | None = None,
    ) -> None:
        super().__init__(name)
        self.x_axis = _vec3(x_axis)
        self.y_axis = None if y_axis is None else _vec3(y_axis)
        self.z_axis = None if z_axis is None else _vec3(z_axis)

    def to_femaster(self) -> str:
        values: list[float] = list(self.x_axis)
        if self.y_axis is not None:
            values.extend(self.y_axis)
        if self.z_axis is not None:
            values.extend(self.z_axis)
        return block([
            keyword("ORIENTATION", NAME=self.name, TYPE="RECTANGULAR"),
            csv(values),
        ])


class CylindricalCoordinateSystem(CoordinateSystem):
    """Cylindrical coordinate system defined by origin, axis and reference."""

    def __init__(
        self,
        name: str,
        origin: Iterable[float],
        axis: Iterable[float],
        reference: Iterable[float],
    ) -> None:
        super().__init__(name)
        self.origin    = _vec3(origin)
        self.axis      = _vec3(axis)
        self.reference = _vec3(reference)

    def to_femaster(self) -> str:
        return block([
            keyword("ORIENTATION", NAME=self.name, TYPE="CYLINDRICAL"),
            csv((*self.origin, *self.axis, *self.reference)),
        ])


class CoordinateSystemRepository(NamedRepository[CoordinateSystem]):
    """Repository of named global coordinate systems."""

    def to_femaster(self) -> str:
        return "\n\n".join(item.to_femaster() for item in self)
