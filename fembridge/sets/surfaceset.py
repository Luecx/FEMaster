from __future__ import annotations
from typing import List, Optional, Iterable, Union

from surfaces.surface import Surface

# Lazy type hints to avoid import cycles
try:
    from elements.element import Element
except Exception:
    Element = object  # type: ignore
try:
    from sets.elementset import ElementSet
except Exception:
    ElementSet = object  # type: ignore

ElemLike = Union[int, "Element"]
ElsetLike = Union[str, "ElementSet"]
SideLike = Union[int, str]


class SurfaceSet:
    """
    Represents a full *SURFACE block (FEMaster), i.e.:

      *SURFACE, TYPE=ELEMENT, NAME=<name>
      <token>, <side>
      <token>, <side>
      ...

    where each line is provided by a `Surface`.

    Only TYPE=ELEMENT is emitted for now (matches your usage).
    """

    def __init__(self, name: str):
        if not name:
            raise ValueError("SurfaceSet name must be non-empty.")
        self.name: str = str(name)
        self._surfaces: List[Surface] = []

    # ------------- building -------------

    def add(self, surface: Surface) -> Surface:
        self._surfaces.append(surface)
        return surface

    def add_element_side(self, elem: ElemLike, side: SideLike) -> Surface:
        s = Surface(elem=elem, side=side)
        self._surfaces.append(s)
        return s

    def add_elset_side(self, elset: ElsetLike, side: SideLike) -> Surface:
        s = Surface(elset=elset, side=side)
        self._surfaces.append(s)
        return s

    def extend(self, items: Iterable[Surface]) -> None:
        for s in items:
            self._surfaces.append(s)

    def __getitem__(self, idx: int) -> Surface:
        return self._surfaces[idx]

    def __len__(self) -> int:
        return len(self._surfaces)

    # ------------- serialization -------------

    def to_femaster(self) -> str:
        """
        Emit:
          *SURFACE, TYPE=ELEMENT, NAME=<name>
          <token>, <side>
          ...
        """
        lines: List[str] = [f"*SURFACE, TYPE=ELEMENT, NAME={self.name}"]
        for s in self._surfaces:
            lines.append(s.to_femaster())
        return "\n".join(lines)

    def to_asami(self) -> str:
        raise NotImplementedError("SurfaceSet.to_asami not implemented yet.")
