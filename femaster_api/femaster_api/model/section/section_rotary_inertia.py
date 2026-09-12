"""Concentrated diagonal rotary inertia on a concrete ``ElementRegion``.

FEMaster's current point-property path accepts the three principal inertia terms;
products of inertia are exported as zero.  The section keeps its assignment
region as an object and turns it into an ``ELSET`` name only when serialized.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..region.region_element import ElementRegion
from .section import Section


class RotaryInertiaSection(Section):
    """Diagonal concentrated rotary inertia property."""

    def __init__(
        self,
        name: str,
        element_region: ElementRegion,
        inertia: Iterable[float],
    ) -> None:
        super().__init__(name, element_region)
        self.inertia = tuple(float(value) for value in inertia)
        if len(self.inertia) != 3:
            raise ValueError(
                "RotaryInertiaSection requires exactly 3 diagonal moments"
            )

    def export(self) -> str:
        return block([
            keyword("ROTARY INERTIA", ELSET=self.element_region.name),
            csv((*self.inertia, 0.0, 0.0, 0.0)),
        ])
