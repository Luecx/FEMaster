"""Named beam profile."""

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject


class Profile(NamedObject):
    """General beam profile using native scalar section properties."""

    def __init__(
        self,
        name: str,
        area: float,
        iy: float,
        iz: float,
        j: float,
        iyz: float = 0.0,
        ey: float = 0.0,
        ez: float = 0.0,
        refy: float = 0.0,
        refz: float = 0.0,
    ) -> None:
        super().__init__(name)
        self.values = tuple(float(value) for value in (
            area, iy, iz, j, iyz, ey, ez, refy, refz
        ))

    def export(self) -> str:
        return block([keyword("PROFILE", NAME=self.name), csv(self.values)])
