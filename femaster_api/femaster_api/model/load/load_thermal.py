"""Structural thermal load referencing a concrete global ``Field`` object.

The thermal load does not duplicate nodal temperatures.  It stores the actual
``Field`` that provides temperature data together with the structural reference
temperature.  The field's semantic name is read only while serializing native
``*TLOAD`` syntax.
"""

from __future__ import annotations

from ..common.format import keyword
from ..field.field import Field
from .load import Load


class ThermalLoad(Load):
    """Structural temperature loading through one concrete field."""

    def __init__(
        self,
        temperature_field: Field,
        reference_temperature: float = 0.0,
    ) -> None:
        if not isinstance(temperature_field, Field):
            raise TypeError("temperature_field must be a Field object")
        self.temperature_field = temperature_field
        self.reference_temperature = float(reference_temperature)

    def export(self, collector: str) -> str:
        return keyword(
            "TLOAD",
            LOAD_COLLECTOR=collector,
            TEMPERATUREFIELD=self.temperature_field.name,
            REFERENCETEMPERATURE=self.reference_temperature,
        )
