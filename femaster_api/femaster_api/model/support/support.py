"""Prescribed structural degrees of freedom on a concrete node-domain target.

A ``Support`` stores a ``Node`` or ``NodeRegion`` object directly together with
up to six prescribed translational/rotational values.  ``None`` leaves a DOF
unconstrained.  An optional orientation is a concrete ``CoordinateSystem``;
neither target IDs/names nor orientation names are accepted as public reference
surrogates.

The owning ``SupportCollector`` supplies only its own name during native export.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..coordinate_system.coordinate_system import CoordinateSystem
from ..node.node import Node
from ..region.region_node import NodeRegion


class Support:
    """One structural support definition owned by a support collector."""

    def __init__(
        self,
        target: Node | NodeRegion,
        values: Iterable[float | None],
        *,
        orientation: CoordinateSystem | None = None,
    ) -> None:
        if not isinstance(target, (Node, NodeRegion)):
            raise TypeError("target must be Node or NodeRegion")
        if orientation is not None and not isinstance(orientation, CoordinateSystem):
            raise TypeError("orientation must be a CoordinateSystem object")

        self.target = target
        self.values = tuple(
            None if value is None else float(value)
            for value in values
        )
        self.orientation = orientation

        if len(self.values) > 6:
            raise ValueError("Support accepts at most 6 prescribed DOF values")

    def export(self, collector: str) -> str:
        """Export this support with the owning collector name."""

        values = list(self.values)
        while values and values[-1] is None:
            values.pop()

        target = self.target.id if isinstance(self.target, Node) else self.target.name
        return block([
            keyword(
                "SUPPORT",
                SUPPORT_COLLECTOR=collector,
                ORIENTATION=(
                    self.orientation.name if self.orientation is not None else None
                ),
            ),
            csv((target, *values)),
        ])
