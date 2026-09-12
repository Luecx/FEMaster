"""Scalar pressure on a concrete surface or surface-region object.

``PressureLoad`` keeps the target and optional amplitude as real model objects.
No public constructor accepts a surface name or amplitude name as a substitute.
The native ``*PLOAD`` record derives those semantic names only at export time.
"""

from __future__ import annotations

from ..amplitude.amplitude import Amplitude
from ..common.format import block, csv, keyword
from ..region.region_surface import SurfaceRegion
from ..surface.surface import Surface
from .load import Load


class PressureLoad(Load):
    """Scalar pressure applied to one concrete surface-domain target."""

    def __init__(
        self,
        target: Surface | SurfaceRegion,
        pressure: float,
        *,
        amplitude: Amplitude | None = None,
    ) -> None:
        if not isinstance(target, (Surface, SurfaceRegion)):
            raise TypeError("target must be Surface or SurfaceRegion")
        if amplitude is not None and not isinstance(amplitude, Amplitude):
            raise TypeError("amplitude must be an Amplitude object")

        self.target = target
        self.pressure = float(pressure)
        self.amplitude = amplitude

    def export(self, collector: str) -> str:
        return block([
            keyword(
                "PLOAD",
                LOAD_COLLECTOR=collector,
                AMPLITUDE=self.amplitude.name if self.amplitude is not None else None,
            ),
            csv((self.target.name, self.pressure)),
        ])
