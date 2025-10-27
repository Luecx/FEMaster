from __future__ import annotations

from typing import Iterable, Optional, Tuple, Union

from ..elements.element import Element
from ..sets.elementset import ElementSet

from .load_base import Load, CoordinateSystem, AmplitudeFunction, TIME_DEPENDENT_ERROR

Target = Union[Element, ElementSet]


class VLoad(Load):
    """Vector load applied to an Element or ElementSet."""

    load_type = "VLOAD"

    def __init__(
        self,
        target: Target,
        vec: Iterable[float],
        coordinate_system: Optional[CoordinateSystem] = None,
        amplitude: Optional[AmplitudeFunction] = None,
    ) -> None:
        super().__init__(coordinate_system=coordinate_system, amplitude=amplitude)

        if not isinstance(target, (Element, ElementSet)):
            raise TypeError("VLoad target must be an Element, ElementSet")

        if isinstance(target, Element):
            self.region = ElementSet.internal(target)
        else:
            self.region = target

        values = tuple(float(x) for x in vec)
        if len(values) != 3:
            raise ValueError("VLoad vector must have length 3 (vx, vy, vz).")

        self.vec: Tuple[float, float, float] = values

    def to_femaster(self, collector_name: str) -> str:
        if self.is_time_dependent:
            raise NotImplementedError(TIME_DEPENDENT_ERROR)
        vx, vy, vz = self.vec
        token = self.region.name if self.region.name else str(int(self.region.elements[0].element_id))
        payload = f"{token}, {vx}, {vy}, {vz}"
        return "\n".join((self._header(collector_name), payload))

    def to_asami(self) -> str:
        raise NotImplementedError("VLoad.to_asami not implemented yet.")
