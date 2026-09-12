"""Base class for one collector-owned FEMaster load definition.

Loads do not live in a second global repository.  Every concrete load instance
belongs directly to one ``LoadCollector`` and receives that collector's semantic
name during export.  This ownership model mirrors how analysis steps activate
groups of loads and avoids duplicated references to the same load object.

Concrete subclasses own their native keyword syntax; the base class defines only
the common export contract.
"""

from __future__ import annotations


class Load:
    """Base class for one load owned by a ``LoadCollector``."""

    def export(self, collector: str) -> str:
        """Export this load using the owning collector name."""

        raise NotImplementedError
