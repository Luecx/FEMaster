"""Concentrated force/moment load on a concrete node-domain target.

``NodalForce`` stores a ``Node`` or ``NodeRegion`` object directly.  Optional
orientation and amplitude relationships likewise store ``CoordinateSystem`` and
``Amplitude`` objects rather than their names.  Native IDs/names are recovered
only in ``export`` when the ``*CLOAD`` record is serialized.

This keeps the editable model free of dangling semantic strings and makes every
legal relationship visible in the public constructor signature.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..amplitude.amplitude import Amplitude
from ..common.format import block, csv, keyword
from ..coordinate_system.coordinate_system import CoordinateSystem
from ..node.node import Node
from ..region.region_node import NodeRegion
from .load import Load


class NodalForce(Load):
    """Six-component concentrated structural load."""

    def __init__(
        self,
        target: Node | NodeRegion,
        values: Iterable[float] = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        *,
        orientation: CoordinateSystem | None = None,
        amplitude: Amplitude | None = None,
    ) -> None:
        if not isinstance(target, (Node, NodeRegion)):
            raise TypeError("target must be Node or NodeRegion")
        if orientation is not None and not isinstance(orientation, CoordinateSystem):
            raise TypeError("orientation must be a CoordinateSystem object")
        if amplitude is not None and not isinstance(amplitude, Amplitude):
            raise TypeError("amplitude must be an Amplitude object")

        self.target = target
        self.values = tuple(float(value) for value in values)
        self.orientation = orientation
        self.amplitude = amplitude

        if len(self.values) != 6:
            raise ValueError(
                "NodalForce requires exactly 6 force/moment values"
            )

    def export(self, collector: str) -> str:
        target = self.target.id if isinstance(self.target, Node) else self.target.name
        return block([
            keyword(
                "CLOAD",
                LOAD_COLLECTOR=collector,
                ORIENTATION=(
                    self.orientation.name if self.orientation is not None else None
                ),
                AMPLITUDE=self.amplitude.name if self.amplitude is not None else None,
            ),
            csv((target, *self.values)),
        ])
