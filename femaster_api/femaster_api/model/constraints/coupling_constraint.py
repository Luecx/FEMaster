"""Coupling constraint data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.sets import EntitySet, NodeSet

from .coupling_type import CouplingType


@dataclass(frozen=True, slots=True)
class CouplingConstraint:
    """Couple a master node set to a slave node, element, or surface set."""

    master: NodeSet
    slave: EntitySet
    dofs: tuple[bool, bool, bool, bool, bool, bool]
    type: CouplingType = CouplingType.KINEMATIC
