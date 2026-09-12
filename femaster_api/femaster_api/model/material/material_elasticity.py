"""Base class for exportable elastic constitutive definitions.

Elasticity is modeled as a material sub-definition rather than as a serializer
case.  Every concrete constitutive law owns the exact native ``*ELASTIC`` block
it requires, while ``Material`` only composes those independent physical
properties in deck order.
"""

from __future__ import annotations


class Elasticity:
    """Base class for one elastic constitutive definition."""

    def export(self) -> str:
        """Export the concrete ``*ELASTIC`` representation."""

        raise NotImplementedError
