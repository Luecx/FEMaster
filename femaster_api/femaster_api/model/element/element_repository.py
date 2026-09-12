"""Repository and grouped export for part-local finite elements.

Element IDs remain sparse semantic FEMaster identifiers.  During export the
repository groups elements by their concrete native ``TYPE`` while preserving
first-seen type order and insertion order inside each group.  This provides a
compact valid deck without moving serializer decisions out of the model layer.
"""

from __future__ import annotations

from ..common.format import block, keyword
from ..common.id_repository import IdRepository
from .element import Element


class ElementRepository(IdRepository[Element]):
    """Own part-local elements keyed by their persistent FEMaster element ID."""

    def export(self) -> str:
        """Export elements as deterministic ``*ELEMENT, TYPE=...`` blocks."""

        groups: dict[str, list[Element]] = {}
        for element in self:
            groups.setdefault(element.type_name, []).append(element)

        rendered: list[str] = []
        for type_name, elements in groups.items():
            rendered.append(
                block([
                    keyword("ELEMENT", TYPE=type_name),
                    *(element.export() for element in elements),
                ])
            )
        return "\n\n".join(rendered)
