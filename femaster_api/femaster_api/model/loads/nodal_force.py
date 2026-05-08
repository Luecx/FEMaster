"""Nodal force model object."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from femaster_api.model.loads.amplitude import Amplitude
from femaster_api.model.nodes import Node
from femaster_api.model.orientations import Orientation
from femaster_api.model.sets import NodeSet

NodalTarget = Node | NodeSet
Vector6 = tuple[float | None, float | None, float | None, float | None, float | None, float | None]


@dataclass(frozen=True, slots=True, init=False)
class NodalForce:
    """Concentrated force or moment applied to a node or node set."""

    target: NodalTarget
    values: Vector6
    orientation: Orientation | None = None
    amplitude: Amplitude | None = None

    def __init__(
        self,
        target: NodalTarget,
        *,
        fx: float | None = None,
        fy: float | None = None,
        fz: float | None = None,
        mx: float | None = None,
        my: float | None = None,
        mz: float | None = None,
        orientation: Orientation | None = None,
        amplitude: Amplitude | None = None,
        values: Iterable[float | None] | None = None,
    ) -> None:
        if values is not None:
            force_values = tuple(None if value is None else float(value) for value in values)
            if len(force_values) != 6:
                raise ValueError("nodal force values must contain six entries")
        else:
            force_values = _six_values(fx, fy, fz, mx, my, mz)
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "values", force_values)
        object.__setattr__(self, "orientation", orientation)
        object.__setattr__(self, "amplitude", amplitude)


def _six_values(*values: float | None) -> Vector6:
    return tuple(None if value is None else float(value) for value in values)  # type: ignore[return-value]
