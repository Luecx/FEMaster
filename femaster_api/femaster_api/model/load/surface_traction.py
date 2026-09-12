"""Distributed surface traction."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.typing import EntityReference
from .load import Load


class SurfaceTraction(Load):
    """Distributed traction vector applied to a surface region."""

    def __init__(
        self,
        target: EntityReference,
        values: Iterable[float],
        *,
        orientation: str | None = None,
        amplitude: str | None = None,
    ) -> None:
        self.target = target
        self.values = tuple(float(value) for value in values)
        self.orientation = orientation
        self.amplitude = amplitude
        if len(self.values) != 3:
            raise ValueError("SurfaceTraction requires exactly 3 components")

    def export(self, collector: str) -> str:
        return block([
            keyword(
                "DLOAD",
                LOAD_COLLECTOR=collector,
                ORIENTATION=self.orientation,
                AMPLITUDE=self.amplitude,
            ),
            csv((self.target, *self.values)),
        ])
