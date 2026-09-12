"""Named repository of section and point-property assignments.

A normal ``SectionRepository`` belongs to a ``Part`` and can contain all section
types.  Every section owns its own native export representation, so the
repository only preserves deterministic order and delegates serialization.
"""

from __future__ import annotations

from ..common.named_repository import NamedRepository
from .section import Section


class SectionRepository(NamedRepository[Section]):
    """Own named section assignments in one part-local scope."""

    def export(self) -> str:
        return "\n\n".join(section.export() for section in self)
