"""Concentrated nodal force and moment load.

``NodalForce`` stores the six structural force/moment components applied to a
node or node region.  Optional coordinate-system and amplitude names remain
semantic references and are emitted directly on the native ``*CLOAD`` keyword.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.typing import EntityReference
from .load import Load


class NodalForce(Load):
    """Six-component concentrated structural load."""

    def __init__(
        self,
        target: EntityReference,
        values: Iterable[float] = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        *,
        orientation: str | None = None,
        amplitude: str | None = None,
    ) -> None:
        self.target = target
        self.values = tuple(float(value) for value in values)
        self.orientation = orientation
        self.amplitude = amplitude
        if len(self.values) != 6:
            raise ValueError(
                "NodalForce requires exactly 6 force/moment values"
            )

    def export(self, collector: str) -> str:
        return block([
            keyword(
                "CLOAD",
                LOAD_COLLECTOR=collector,
                ORIENTATION=self.orientation,
                AMPLITUDE=self.amplitude,
            ),
            csv((self.target, *self.values)),
        ])
