"""Support model object."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from femaster_api.model.nodes import Node
from femaster_api.model.orientations import Orientation
from femaster_api.model.sets import NodeSet

NodalTarget = Node | NodeSet
SupportValues = tuple[float | None, float | None, float | None, float | None, float | None, float | None]


@dataclass(frozen=True, slots=True, init=False)
class Support:
    """Prescribed nodal degrees of freedom for a node or node set."""

    target: NodalTarget
    values: SupportValues
    orientation: Orientation | None = None

    def __init__(
        self,
        target: NodalTarget,
        *,
        ux: float | None = None,
        uy: float | None = None,
        uz: float | None = None,
        rx: float | None = None,
        ry: float | None = None,
        rz: float | None = None,
        orientation: Orientation | None = None,
        values: Iterable[float | None] | None = None,
    ) -> None:
        if values is not None:
            support_values = tuple(None if value is None else float(value) for value in values)
            if len(support_values) != 6:
                raise ValueError("support values must contain six entries")
        else:
            support_values = _six_values(ux, uy, uz, rx, ry, rz)
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "values", support_values)
        object.__setattr__(self, "orientation", orientation)


def _six_values(*values: float | None) -> SupportValues:
    return tuple(None if value is None else float(value) for value in values)  # type: ignore[return-value]
