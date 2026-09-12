"""Element-region specialization storing actual finite-element objects.

``ElementRegion`` corresponds to native ``*ELSET`` and owns references to
``Element`` instances rather than sparse integer IDs.  The persistent element
IDs remain part of the element objects and are emitted only during serialization.
Sections, loads and constraints can consequently hold the region object itself.
"""

from __future__ import annotations

from ..element.element import Element
from .region import Region


class ElementRegion(Region[Element]):
    """Ordered region of concrete elements exported through ``*ELSET``."""

    keyword_name = "ELSET"
    name_key = "NAME"
    member_type = Element

    def _member_value(self, member: Element) -> int:
        return member.id
