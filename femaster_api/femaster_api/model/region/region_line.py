"""Concrete line region.

Semantic region of line entities; FEMaster currently has no standalone line-set keyword.  The class contains only domain-specific keyword metadata; member
storage, deterministic ordering and row export are implemented by ``Region``.
Keeping each region domain in its own module makes imports and future
domain-specific validation explicit.
"""

from __future__ import annotations

from .region import Region


class LineRegion(Region):
    """Semantic region of line entities; FEMaster currently has no standalone line-set keyword."""

    keyword_name = None
    name_key = 'NAME'
