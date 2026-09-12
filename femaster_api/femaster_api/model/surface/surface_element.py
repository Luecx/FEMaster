"""Element-side surface definition using concrete element-model objects.

``ElementSurface`` represents native ``*SURFACE, TYPE=ELEMENT``.  Each entry is
an ``Element`` or ``ElementRegion`` together with the addressed boundary side;
integer IDs and region names are serialization details only.  This makes surface
topology directly navigable and prevents dangling text references in the Python
model.

The class remains distinct from ``SurfaceRegion``: a surface defines topology,
whereas a surface region groups already-defined surfaces.
"""

from __future__ import annotations

from ..common.format import block, csv, keyword
from ..element.element import Element
from ..region.region_element import ElementRegion
from .surface import Surface


class ElementSurface(Surface):
    """Named surface assembled from elements or element regions."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.entries: list[tuple[Element | ElementRegion, int]] = []

    def add(
        self,
        element: Element | ElementRegion,
        side: int,
    ) -> "ElementSurface":
        """Append one concrete element-side target and return this surface."""

        if not isinstance(element, (Element, ElementRegion)):
            raise TypeError("surface target must be Element or ElementRegion")

        side = int(side)
        if side <= 0:
            raise ValueError("surface side must be positive")

        self.entries.append((element, side))
        return self

    def export(self) -> str:
        """Export one ``*SURFACE, TYPE=ELEMENT`` block."""

        def value(target: Element | ElementRegion) -> int | str:
            return target.id if isinstance(target, Element) else target.name

        return block([
            keyword("SURFACE", NAME=self.name, TYPE="ELEMENT"),
            *(csv((value(target), f"S{side}")) for target, side in self.entries),
        ])
