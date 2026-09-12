"""Rotary-inertia point property."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .section import Section


class RotaryInertiaSection(Section):
    """Diagonal concentrated rotary inertia assigned to ROTARYI elements."""

    def __init__(self, name: str, element_region: str, inertia: Iterable[float]) -> None:
        super().__init__(name, element_region)
        self.inertia = tuple(float(value) for value in inertia)
        if len(self.inertia) != 3:
            raise ValueError("RotaryInertiaSection requires exactly 3 diagonal moments")

    def export(self) -> str:
        return block([
            keyword("ROTARY INERTIA", ELSET=self.element_region),
            csv((*self.inertia, 0.0, 0.0, 0.0)),
        ])
