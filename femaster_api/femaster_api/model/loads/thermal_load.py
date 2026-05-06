"""Thermal load data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.fields import Field


@dataclass(frozen=True, slots=True)
class ThermalLoad:
    temperature_field: Field
    reference_temperature: float
