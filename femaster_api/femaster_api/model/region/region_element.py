"""Concrete element region.

Region of elements exported through the native ``*ELSET`` keyword.  The class contains only domain-specific keyword metadata; member
storage, deterministic ordering and row export are implemented by ``Region``.
Keeping each region domain in its own module makes imports and future
domain-specific validation explicit.
"""

from __future__ import annotations

from .region import Region


class ElementRegion(Region):
    """Region of elements exported through the native ``*ELSET`` keyword."""

    keyword_name = 'ELSET'
    name_key = 'NAME'
