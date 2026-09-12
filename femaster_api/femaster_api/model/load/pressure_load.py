"""Scalar surface pressure."""

from ..common.format import block, csv, keyword
from ..common.typing import EntityReference
from .load import Load


class PressureLoad(Load):
    """Scalar pressure applied to a surface region."""

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
            keyword("PLOAD", LOAD_COLLECTOR=collector, AMPLITUDE=self.amplitude),
            csv((self.target, self.pressure)),
        ])
