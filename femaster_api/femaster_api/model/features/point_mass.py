"""Point mass feature data object."""

from __future__ import annotations

from dataclasses import dataclass

from femaster_api.model.sets import NodeSet

Vector3 = tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class PointMass:
    """Lumped mass and optional springs applied to a node set."""

    node_set: NodeSet
    mass: float = 0.0
    inertia: Vector3 = (0.0, 0.0, 0.0)
    spring: Vector3 = (0.0, 0.0, 0.0)
    rotational_spring: Vector3 = (0.0, 0.0, 0.0)
