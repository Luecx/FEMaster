"""Surface-region specialization grouping concrete ``Surface`` definitions.

A ``SurfaceRegion`` maps to native ``*SFSET``.  Its members are the already
materialized surface objects themselves, not their names.  This keeps topology
ownership explicit while still exporting exactly the semantic names expected by
FEMaster's input format.
"""

from __future__ import annotations

from ..surface.surface import Surface
from .region import Region


class SurfaceRegion(Region[Surface]):
    """Ordered group of concrete surfaces exported through ``*SFSET``."""

    keyword_name = "SFSET"
    name_key = "SFSET"
    member_type = Surface

    def _member_value(self, member: Surface) -> str:
        return member.name
