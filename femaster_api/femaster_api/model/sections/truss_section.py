"""Truss section data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.materials import Material
from femaster_api.model.sets import EntitySet


@dataclass(frozen=True, slots=True)
class TrussSection:
    name: str
    material: Material
    element_set: EntitySet
    area: float
