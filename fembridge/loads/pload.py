from __future__ import annotations

from typing import Optional, Union

from ..sets.surfaceset import SurfaceSet

from .load_base import Load, AmplitudeFunction, TIME_DEPENDENT_ERROR


class PLoad(Load):
    """Pressure load applied to a SurfaceSet."""

    load_type = "PLOAD"

    def __init__(
        self,
        surface_set: Union[SurfaceSet, str],
        p: float,
        amplitude: Optional[AmplitudeFunction] = None,
    ) -> None:
        super().__init__(amplitude=amplitude)
        self.name = self._resolve_name(surface_set)
        self.p = float(p)

    @staticmethod
    def _resolve_name(surface_set: Union[SurfaceSet, str]) -> str:
        if isinstance(surface_set, SurfaceSet):
            name = surface_set.name
        else:
            name = str(surface_set).strip()
        if not name:
            raise ValueError("PLoad requires a non-empty SurfaceSet name.")
        return name

    def to_femaster(self, collector_name: str) -> str:
        if self.is_time_dependent:
            raise NotImplementedError(TIME_DEPENDENT_ERROR)
        payload = f"{self.name}, {self.p}"
        return "\n".join((self._header(collector_name), payload))

    def to_asami(self) -> str:
        raise NotImplementedError("PLoad.to_asami not implemented yet.")
