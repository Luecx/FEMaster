from __future__ import annotations

from typing import Iterable, Optional, Tuple

from ..sets.surfaceset import SurfaceSet

from .load_base import Load, CoordinateSystem, AmplitudeFunction, TIME_DEPENDENT_ERROR


class DLoad(Load):
    """Distributed vector load applied to a SurfaceSet."""

    load_type = "DLOAD"

    def __init__(
        self,
        surface_set: SurfaceSet,
        vec: Iterable[float],
        coordinate_system: Optional[CoordinateSystem] = None,
        amplitude: Optional[AmplitudeFunction] = None,
    ) -> None:
        super().__init__(coordinate_system=coordinate_system, amplitude=amplitude)

        if not isinstance(surface_set, SurfaceSet):
            raise TypeError("DLoad surface_set must be a SurfaceSet instance.")
        if not surface_set.name:
            raise ValueError("DLoad requires a SurfaceSet with a non-empty name.")

        self.surface_set = surface_set

        values = tuple(float(x) for x in vec)
        if len(values) != 3:
            raise ValueError("DLoad vector must have length 3 (tx, ty, tz).")
        self.vec: Tuple[float, float, float] = values

    def to_femaster(self, collector_name: str) -> str:
        if self.is_time_dependent:
            raise NotImplementedError(TIME_DEPENDENT_ERROR)
        tx, ty, tz = self.vec
        payload = f"{self.surface_set.name}, {tx}, {ty}, {tz}"
        return "\n".join((self._header(collector_name), payload))

    def to_asami(self) -> str:
        raise NotImplementedError("DLoad.to_asami not implemented yet.")
