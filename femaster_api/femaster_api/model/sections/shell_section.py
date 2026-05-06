"""Shell section data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.materials import Material
from femaster_api.model.orientations import Orientation
from femaster_api.model.sets import EntitySet


@dataclass(frozen=True, slots=True)
class ShellSection:
    name: str
    material: Material
    element_set: EntitySet
    thickness: float
    orientation: Orientation | None = None
