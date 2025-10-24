from __future__ import annotations

from typing import Iterable, Optional, Tuple, Union

from sets.surfaceset import SurfaceSet

from .load_base import Load, CoordinateSystem


class DLoad(Load):
    """Distributed vector load applied to a SurfaceSet."""

    load_type = "DLOAD"

    def __init__(
        self,
        surface_set: Union[SurfaceSet, str],
        vec: Iterable[float],
        coordinate_system: Optional[CoordinateSystem] = None,
    ) -> None:
        super().__init__(coordinate_system)

        self.name = self._resolve_name(surface_set)

        v = tuple(float(x) for x in vec)
        if len(v) != 3:
            raise ValueError("DLoad vector must have length 3 (tx, ty, tz).")
        self.vec: Tuple[float, float, float] = v

    @staticmethod
    def _resolve_name(surface_set: Union[SurfaceSet, str]) -> str:
        if isinstance(surface_set, SurfaceSet):
            name = surface_set.name
        else:
            name = str(surface_set).strip()
        if not name:
            raise ValueError("DLoad requires a non-empty SurfaceSet name.")
        return name

    def to_femaster(self, collector_name: str) -> str:
        tx, ty, tz = self.vec
        payload = f"{self.name}, {tx}, {ty}, {tz}"
        return "\n".join((self._header(collector_name), payload))

    def to_asami(self) -> str:
        raise NotImplementedError("DLoad.to_asami not implemented yet.")
