"""Scalar pressure load on a surface region.

Pressure is modeled separately from general vector traction because the native
``*PLOAD`` syntax and physical interpretation are distinct.  An optional
amplitude reference allows the same pressure definition to participate in
time-dependent analyses.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from ..common.typing import EntityReference
from .load import Load


class PressureLoad(Load):
    """Scalar pressure applied to one surface target."""

    def __init__(
        self,
        target: EntityReference,
        pressure: float,
        *,
        amplitude: str | None = None,
    ) -> None:
        self.target = target
        self.pressure = float(pressure)
        self.amplitude = amplitude

    def export(self, collector: str) -> str:
        return block([
            keyword(
                "PLOAD",
                LOAD_COLLECTOR=collector,
                AMPLITUDE=self.amplitude,
            ),
            csv((self.target, self.pressure)),
        ])
