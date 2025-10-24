from __future__ import annotations

from typing import Iterable, Optional, Tuple, Union

from elements.element import Element
from sets.elementset import ElementSet

from .load_base import Load, CoordinateSystem

Target = Union[Element, ElementSet, int]


class VLoad(Load):
    """Vector load applied to an Element or ElementSet."""

    load_type = "VLOAD"

    def __init__(
        self,
        target: Target,
        vec: Iterable[float],
        coordinate_system: Optional[CoordinateSystem] = None,
    ) -> None:
        super().__init__(coordinate_system)

        if not isinstance(target, (Element, ElementSet, int)):
            raise TypeError("VLoad target must be an Element, ElementSet, or element id (int).")

        self.token = self._resolve_target(target)

        v = tuple(float(x) for x in vec)
        if len(v) != 3:
            raise ValueError("VLoad vector must have length 3 (vx, vy, vz).")
        self.vec: Tuple[float, float, float] = v

    @staticmethod
    def _resolve_target(target: Target) -> str:
        if isinstance(target, Element):
            if target.element_id is None:
                raise ValueError("Element used in VLoad has no element_id.")
            return str(int(target.element_id))
        if isinstance(target, ElementSet):
            if not target.name:
                raise ValueError("ElementSet used in VLoad has empty name.")
            return target.name
        return str(int(target))

    def to_femaster(self, collector_name: str) -> str:
        vx, vy, vz = self.vec
        payload = f"{self.token}, {vx}, {vy}, {vz}"
        return "\n".join((self._header(collector_name), payload))

    def to_asami(self) -> str:
        raise NotImplementedError("VLoad.to_asami not implemented yet.")
