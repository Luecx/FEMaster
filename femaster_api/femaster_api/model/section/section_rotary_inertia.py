"""Concentrated rotary inertia assigned to ``ROTARYI`` point elements.

The native command currently supports the three principal inertia values while
products of inertia are required to remain zero.  The Python model therefore
stores the physically active diagonal components explicitly and emits the three
zero products in FEMaster's six-value input order.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .section import Section


class RotaryInertiaSection(Section):
    """Principal concentrated rotary-inertia property."""

    def __init__(
        self,
        name: str,
        element_region: str,
        inertia: Iterable[float],
    ) -> None:
        super().__init__(name, element_region)
        self.inertia = tuple(float(value) for value in inertia)
        if len(self.inertia) != 3:
            raise ValueError(
                "RotaryInertiaSection requires exactly 3 principal inertias"
            )

    def export(self) -> str:
        return block([
            keyword("ROTARY INERTIA", ELSET=self.element_region),
            csv((*self.inertia, 0.0, 0.0, 0.0)),
        ])
