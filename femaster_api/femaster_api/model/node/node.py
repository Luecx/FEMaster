"""Part-local finite-element node definition.

A ``Node`` stores the sparse identifier and reference coordinates exactly as
supplied by the user or imported deck.  IDs are never replaced by repository
positions.  The owning ``Part`` controls scope; this class only represents one
row of a native ``*NODE`` block and can therefore export its own data row.

Nodes remain intentionally lightweight because all set membership and assembly
placement belong to regions and instances rather than to the node object itself.
"""

from __future__ import annotations

from ..common.format import csv


class Node:
    """One part-local FEMaster node with a persistent integer identifier."""

    def __init__(self, id: int, x: float, y: float, z: float = 0.0) -> None:
        self.id = int(id)
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def export(self) -> str:
        """Export this node as one ``*NODE`` data row."""

        return csv((self.id, self.x, self.y, self.z))
