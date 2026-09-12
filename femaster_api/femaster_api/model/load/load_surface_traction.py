"""Distributed vector traction on a concrete surface-domain target.

The target is a ``Surface`` or ``SurfaceRegion`` object.  Optional orientation
and amplitude relationships are concrete ``CoordinateSystem`` / ``Amplitude``
objects as well.  The native ``*DLOAD`` names are therefore serialization output,
not the in-memory representation of model relationships.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..amplitude.amplitude import Amplitude
from ..common.format import block, csv, keyword
from ..coordinate_system.coordinate_system import CoordinateSystem
from ..region.region_surface import SurfaceRegion
from ..surface.surface import Surface
from .load import Load


class SurfaceTraction(Load):
    """Three-component distributed traction on a surface target."""

    def __init__(
        self,
        target: Surface | SurfaceRegion,
        values: Iterable[float],
        *,
        orientation: CoordinateSystem | None = None,
        amplitude: Amplitude | None = None,
    ) -> None:
        if not isinstance(target, (Surface, SurfaceRegion)):
            raise TypeError("target must be Surface or SurfaceRegion")
        if orientation is not None and not isinstance(orientation, CoordinateSystem):
            raise TypeError("orientation must be a CoordinateSystem object")
        if amplitude is not None and not isinstance(amplitude, Amplitude):
            raise TypeError("amplitude must be an Amplitude object")

        self.target = target
        self.values = tuple(float(value) for value in values)
        self.orientation = orientation
        self.amplitude = amplitude

        if len(self.values) != 3:
            raise ValueError("SurfaceTraction requires exactly 3 components")

    def export(self, collector: str) -> str:
        return block([
            keyword(
                "DLOAD",
                LOAD_COLLECTOR=collector,
                ORIENTATION=(
                    self.orientation.name if self.orientation is not None else None
                ),
                AMPLITUDE=self.amplitude.name if self.amplitude is not None else None,
            ),
            csv((self.target.name, *self.values)),
        ])
