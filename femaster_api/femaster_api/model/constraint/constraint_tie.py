"""Tie relation between a master and slave region.

The object preserves native adjustment and optional search-distance controls.
Master and slave are semantic names because region resolution is performed by
FEMaster after the complete assembly has been built.
"""

from __future__ import annotations

from ..common.format import keyword
from .constraint import Constraint


class Tie(Constraint):
    """Tie one slave region to one master region."""

    def __init__(
        self,
        master: str,
        slave: str,
        *,
        adjust: bool = True,
        distance: float | None = None,
    ) -> None:
        self.master = str(master)
        self.slave = str(slave)
        self.adjust = bool(adjust)
        self.distance = None if distance is None else float(distance)

    def export(self) -> str:
        return keyword(
            "TIE",
            MASTER=self.master,
            SLAVE=self.slave,
            ADJUST="YES" if self.adjust else "NO",
            DISTANCE=self.distance,
        )
