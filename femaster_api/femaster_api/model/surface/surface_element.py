"""Element-side surface definition.

``ElementSurface`` represents the native ``*SURFACE, TYPE=ELEMENT`` form.  Each
entry addresses an element or element region and one boundary side.  The entry
keeps the semantic reference exactly as supplied, so assembly-qualified or named
references are not converted to repository positions.

This class is separate from ``SurfaceRegion``: the surface defines topology,
whereas a surface region groups already-defined surfaces.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from ..common.typing import ElementReference
from .surface import Surface


class ElementSurface(Surface):
    """Named surface assembled from element-side references."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.entries: list[tuple[ElementReference, int]] = []

    def add(self, element: ElementReference, side: int) -> "ElementSurface":
        """Append one element-side reference and return this surface."""

        side = int(side)
        if side <= 0:
            raise ValueError("surface side must be positive")
        self.entries.append((element, side))
        return self

    def export(self) -> str:
        """Export one ``*SURFACE, TYPE=ELEMENT`` block."""

        return block([
            keyword("SURFACE", NAME=self.name, TYPE="ELEMENT"),
            *(csv((element, f"S{side}")) for element, side in self.entries),
        ])
