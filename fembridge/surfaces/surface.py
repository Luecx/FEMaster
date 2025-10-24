from __future__ import annotations
from typing import Union

from elements.element import Element
from sets.elementset import ElementSet

SideLike = Union[int, str]


class Surface:
    """
    Eine einzelne Zeile in einem *SURFACE-Block.

    - `elems`: Element oder ElementSet.
      Intern wird immer ein ElementSet gehalten.
      * Bei Element → temporäres, unbenanntes Set mit genau diesem Element.
      * Bei ElementSet → direkt verwendet.
    - `side`: Seitenangabe (z. B. 1, 2, "SNEG", ...).
    """

    def __init__(self, elems: Union[Element, ElementSet], side: SideLike) -> None:
        if isinstance(elems, ElementSet):
            self._elset: ElementSet = elems
        elif isinstance(elems, Element):
            self._elset = ElementSet.internal(elems)
        else:
            raise TypeError("Surface: 'elems' must be Element or ElementSet.")

        self._side: str = side if isinstance(side, str) else f"S{side}"

    @property
    def elset(self) -> ElementSet:
        return self._elset

    @property
    def side(self) -> str:
        return self._side

    def to_femaster(self) -> str:
        token = self._elset.name or str(self._elset.elements[0].element_id)
        return f"{token}, {self._side}"

    def to_asami(self) -> str:
        raise NotImplementedError("Surface.to_asami not implemented yet.")
