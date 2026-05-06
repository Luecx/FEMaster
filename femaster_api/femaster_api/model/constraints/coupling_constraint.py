"""Coupling constraint data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.sets import EntitySet

from .coupling_type import CouplingType


@dataclass(frozen=True, slots=True)
class CouplingConstraint:
    master: EntitySet
    slave: EntitySet
    dofs: tuple[bool, bool, bool, bool, bool, bool]
    type: CouplingType = CouplingType.KINEMATIC
