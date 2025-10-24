from __future__ import annotations
from typing import Iterable, List, Union, Iterator
from elements.element import Element

ElemListLike = Union[Element, Iterable[Element], "ElementSet"]


class ElementSet:
    """
    Container ausschließlich für Element-Objekte (niemals IDs).
    name kann None sein (interne/temporäre Sets).
    """

    def __init__(self, name: str | None = None, elements: List[Element] | None = None) -> None:
        self.name: str | None = (str(name) if name else None)
        self.elements: List[Element] = []
        if elements:
            for e in elements:
                if not isinstance(e, Element):
                    raise TypeError("ElementSet darf nur Element-Instanzen enthalten.")
            self.elements.extend(elements)

    # ---------------- Fabriken ----------------

    @staticmethod
    def internal(value: ElemListLike) -> ElementSet:
        """Unbenanntes Set aus Element | Iterable[Element] | ElementSet."""
        if isinstance(value, ElementSet):
            return value
        return ElementSet(None, ElementSet._normalize(value))

    @staticmethod
    def create(name: str, value: ElemListLike) -> ElementSet:
        """Benanntes Set (kopiert Elemente falls ElementSet)."""
        if not name:
            raise ValueError("ElementSet.create(): name darf nicht leer sein.")
        if isinstance(value, ElementSet):
            return ElementSet(name, list(value.elements))
        return ElementSet(name, ElementSet._normalize(value))

    # --------------- Normalisierung ---------------

    @staticmethod
    def _normalize(value: ElemListLike) -> List[Element]:
        """Konvertiert Element | Iterable[Element] | ElementSet -> list[Element]."""
        if isinstance(value, Element):
            return [value]
        if isinstance(value, ElementSet):
            return list(value.elements)

        try:
            lst = list(value)  # type: ignore[arg-type]
        except TypeError:
            raise TypeError("Erwarte Element, Iterable[Element] oder ElementSet.")

        for e in lst:
            if not isinstance(e, Element):
                raise TypeError("Iterable muss ausschließlich Element-Instanzen enthalten.")
        return lst

    # ---------------- Convenience ----------------

    def add(self, value: ElemListLike) -> None:
        """Fügt Element | Iterable[Element] | ElementSet hinzu."""
        self.elements.extend(ElementSet._normalize(value))

    def __iter__(self) -> Iterator[Element]:
        return iter(self.elements)

    def __len__(self) -> int:
        return len(self.elements)

    # ---------------- Exports ----------------

    def to_femaster(self) -> str:
        if not self.name:
            raise ValueError("Unbenanntes ElementSet kann nicht exportiert werden. Verwende ElementSet.create(name, ...).")

        ids: List[int] = []
        for e in self.elements:
            eid = getattr(e, "element_id", None)
            if eid is None:
                raise ValueError(f"ElementSet '{self.name}' enthält Element ohne ID.")
            ids.append(int(eid))

        lines = [f"*ELSET, NAME={self.name}"]
        line: List[str] = []
        for i, eid in enumerate(ids, 1):
            line.append(str(eid))
            if i % 16 == 0:
                lines.append(", ".join(line))
                line = []
        if line:
            lines.append(", ".join(line))
        return "\n".join(lines)

    def to_asami(self) -> str:
        raise NotImplementedError
