"""Tie relation between concrete master and slave surface-domain objects.

``Tie`` stores ``Surface`` or ``SurfaceRegion`` objects for both sides.  Native
surface names are generated only during export, so the constraint cannot retain
unresolved master/slave strings.  Adjustment and optional search distance remain
intrinsic tie parameters.
"""

from __future__ import annotations

from ..common.format import keyword
from ..region.region_surface import SurfaceRegion
from ..surface.surface import Surface
from .constraint import Constraint


class Tie(Constraint):
    """Tie one slave surface target to one master surface target."""

    def __init__(
        self,
        master: Surface | SurfaceRegion,
        slave: Surface | SurfaceRegion,
        *,
        adjust: bool = True,
        distance: float | None = None,
    ) -> None:
        if not isinstance(master, (Surface, SurfaceRegion)):
            raise TypeError("master must be Surface or SurfaceRegion")
        if not isinstance(slave, (Surface, SurfaceRegion)):
            raise TypeError("slave must be Surface or SurfaceRegion")

        self.master = master
        self.slave = slave
        self.adjust = bool(adjust)
        self.distance = None if distance is None else float(distance)

    def export(self) -> str:
        return keyword(
            "TIE",
            MASTER=self.master.name,
            SLAVE=self.slave.name,
            ADJUST="YES" if self.adjust else "NO",
            DISTANCE=self.distance,
        )
