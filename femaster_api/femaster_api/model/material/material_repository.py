"""Named repository of project-level materials.

Materials are globally shared definitions, so one ``Project`` owns one
``MaterialRepository``.  The repository adds no constitutive logic; it only
provides deterministic name/position lookup and delegates native export to each
material object.
"""

from __future__ import annotations

from ..common.named_repository import NamedRepository
from .material import Material


class MaterialRepository(NamedRepository[Material]):
    """Own globally named materials in deterministic insertion order."""

    def export(self) -> str:
        return "\n\n".join(material.export() for material in self)
