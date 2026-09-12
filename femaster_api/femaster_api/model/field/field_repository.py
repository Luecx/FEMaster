"""Named repository of globally defined model fields.

Project-level fields may be referenced by thermal loads, topology workflows or
other analyses.  Result frames use their own direct field dictionary because
result field names are scoped by frame rather than by the editable project.
"""

from __future__ import annotations

from ..common.named_repository import NamedRepository
from .field import Field


class FieldRepository(NamedRepository[Field]):
    """Own named model fields and export them in deterministic order."""

    def export(self) -> str:
        return "\n\n".join(field.export() for field in self)
