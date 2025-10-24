from __future__ import annotations

from typing import Iterable, Tuple

from .coordinate_system_base import CoordinateSystem


class CylindricalCoordinateSystem(CoordinateSystem):
    """Cylindrical coordinate system defined by base, radial, and theta points."""

    def __init__(
        self,
        name: str,
        base_point: Iterable[float],
        radial_point: Iterable[float],
        theta_point: Iterable[float],
    ) -> None:
        super().__init__(name, "CYLINDRICAL", base_point)
        self.radial_point: Tuple[float, float, float] = self._to_triplet(radial_point, "radial_point")
        self.theta_point: Tuple[float, float, float] = self._to_triplet(theta_point, "theta_point")

    def to_femaster(self) -> str:
        lines = [
            self._header(),
            self._format_triplet(self.base_point),
            self._format_triplet(self.radial_point),
            self._format_triplet(self.theta_point),
        ]
        return "\n".join(lines)
