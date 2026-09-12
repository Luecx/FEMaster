"""Named repository of global coordinate systems.

The repository exists at project scope and delegates all orientation-specific
syntax to the concrete coordinate-system classes.  Names remain the persistent
identity used by sections, loads and supports.
"""

from __future__ import annotations

from ..common.named_repository import NamedRepository
from .coordinate_system import CoordinateSystem


class CoordinateSystemRepository(NamedRepository[CoordinateSystem]):
    """Own globally named coordinate systems."""

    def export(self) -> str:
        return "\n\n".join(item.export() for item in self)
