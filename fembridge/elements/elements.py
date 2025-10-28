from __future__ import annotations
from functools import lru_cache
from typing import Dict, List, Optional
from .element import Element

class Elements:
    """
    Sparse 0-based list of Element objects.
    IDs are assigned automatically if Element.element_id is None.
    FEMaster output is grouped by element type:
      *ELEMENT, TYPE=XXXX
      id, n1, n2, ...
    """
    def __init__(self) -> None:
        self._items: List[Optional[Element]] = []

    def add(self, elem: Element) -> Element:
        if elem.element_id is None:
            for i, v in enumerate(self._items):
                if v is None:
                    elem.element_id = i
                    self._items[i] = elem
                    return elem
            elem.element_id = len(self._items)
            self._items.append(elem)
            return elem
        else:
            eid = int(elem.element_id)
            if eid < 0:
                raise ValueError("element_id must be >= 0")
            if eid >= len(self._items):
                self._items.extend([None] * (eid + 1 - len(self._items)))
            self._items[eid] = elem
            return elem

    def __getitem__(self, element_id: int) -> Element:
        val = self._items[element_id]
        if val is None:
            raise KeyError(f"Element id {element_id} is empty.")
        return val

    def __len__(self) -> int:
        return len(self._items)

    def to_femaster(self) -> str:
        # Preserve first-appearance order of types
        order: List[str] = []
        buckets: Dict[str, List[Element]] = {}
        for e in self._items:
            if e is None:
                continue
            t = e.elem_type
            if t not in buckets:
                buckets[t] = []
                order.append(t)
            buckets[t].append(e)

        lines: List[str] = []
        for t in order:
            lines.append(f"*ELEMENT, TYPE={t}")
            for e in buckets[t]:
                lines.append(e.femaster_record())
        return "\n".join(lines)

    def to_asami(self) -> str:
        raise NotImplementedError("Elements.to_asami not implemented yet.")

    # --- helpers for geometry utilities -----------------------------------
    @staticmethod
    @lru_cache(maxsize=None)
    def _element_type_map() -> Dict[str, type[Element]]:
        from . import (
            C2D3,
            C2D4,
            C2D6,
            C2D8,
            S3,
            S4,
            S6,
            S8,
            C3D4,
            C3D5,
            C3D6,
            C3D8,
            C3D10,
            C3D15,
            C3D20,
        )

        return {
            "C2D3": C2D3,
            "C2D4": C2D4,
            "C2D6": C2D6,
            "C2D8": C2D8,
            "S3": S3,
            "S4": S4,
            "S6": S6,
            "S8": S8,
            "C3D4": C3D4,
            "C3D5": C3D5,
            "C3D6": C3D6,
            "C3D8": C3D8,
            "C3D10": C3D10,
            "C3D15": C3D15,
            "C3D20": C3D20,
        }

    @classmethod
    def element_class_for(cls, elem_type: str) -> type[Element]:
        try:
            return cls._element_type_map()[elem_type.upper()]
        except KeyError as exc:
            raise ValueError(f"Unsupported element type '{elem_type}'.") from exc

    def determine_element_order(self) -> int:
        first_order = {
            "C2D3",
            "C2D4",
            "C3D4",
            "C3D6",
            "C3D8",
            "S3",
            "S4",
            "S6",
            "S8",
        }
        second_order = {"C2D6", "C2D8", "C3D10", "C3D15", "C3D20"}

        for element in self._items:
            if element is None:
                continue
            etype = element.elem_type
            if etype in first_order:
                return 1
            if etype in second_order:
                return 2
        return 0
