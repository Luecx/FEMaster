"""Tie constraint."""

from ..common.format import keyword
from .constraint import Constraint


class Tie(Constraint):
    """Tie slave nodes/surfaces to a master surface or line."""

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
