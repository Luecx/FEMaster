"""Concrete surface region.

Region of named/materialized surfaces exported through ``*SFSET``.  The class contains only domain-specific keyword metadata; member
storage, deterministic ordering and row export are implemented by ``Region``.
Keeping each region domain in its own module makes imports and future
domain-specific validation explicit.
"""

from __future__ import annotations

from .region import Region


class SurfaceRegion(Region):
    """Region of named/materialized surfaces exported through ``*SFSET``."""

    keyword_name = 'SFSET'
    name_key = 'SFSET'
