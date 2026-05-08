"""Beam section data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.materials import Material
from femaster_api.model.sets import ElementSet

from .profile import Profile

Vector3 = tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class BeamSection:
    """Beam section assignment for an element set."""

    name: str
    material: Material
    element_set: ElementSet
    profile: Profile
    orientation: Vector3 | None = None
