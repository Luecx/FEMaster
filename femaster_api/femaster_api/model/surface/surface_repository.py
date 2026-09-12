"""Named repository for surface definitions in one FEMaster scope.

A ``Part`` and the assembled ``Project`` each own an independent
``SurfaceRepository``.  The repository provides name/position lookup and writes
all contained surfaces in deterministic insertion order.  Concrete surface
syntax remains on the surface class itself.
"""

from __future__ import annotations

from ..common.named_repository import NamedRepository
from .surface import Surface


class SurfaceRepository(NamedRepository[Surface]):
    """Own named surfaces and delegate export to each concrete surface."""

    def export(self) -> str:
        """Export all surfaces in deterministic repository order."""

        return "\n\n".join(surface.export() for surface in self)
