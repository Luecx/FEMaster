"""Structural thermal load referencing a named temperature field.

The thermal load does not duplicate nodal temperatures.  It references a global
``Field`` by name and stores only the structural reference temperature used to
convert temperature differences into thermal strain.
"""

from __future__ import annotations

from ..common.format import keyword
from .load import Load


class ThermalLoad(Load):
    """Structural temperature loading through a named temperature field."""

    def __init__(
        self,
        temperature_field: str,
        reference_temperature: float = 0.0,
    ) -> None:
        self.temperature_field = str(temperature_field)
        self.reference_temperature = float(reference_temperature)

    def export(self, collector: str) -> str:
        return keyword(
            "TLOAD",
            LOAD_COLLECTOR=collector,
            TEMPERATUREFIELD=self.temperature_field,
            REFERENCETEMPERATURE=self.reference_temperature,
        )
