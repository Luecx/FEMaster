"""Semantic line-element region without a standalone native set keyword.

FEMaster currently exposes no independent line-entity object or native line-set
keyword.  The Python representation therefore groups concrete ``Element``
objects that are interpreted as line-like by the consuming operation.  Keeping
real element objects still follows the same no-string-reference rule as every
other modeled relationship.
"""

from __future__ import annotations

from ..element.element import Element
from .region import Region


class LineRegion(Region[Element]):
    """Semantic region of line-like finite elements."""

    keyword_name = None
    name_key = "NAME"
    member_type = Element

    def _member_value(self, member: Element) -> int:
        return member.id
