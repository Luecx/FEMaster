from __future__ import annotations

from typing import Iterable, Optional, Tuple

from .coordinate_system_base import CoordinateSystem


class RectangularCoordinateSystem(CoordinateSystem):
    """
    Rectangular coordinate system definition.

    Supports providing between one and three direction vectors. Any additional
    vectors beyond the first are simply stored and emitted if supplied, but not
    required by FEMaster.
    """

    def __init__(
        self,
        name: str,
        base_point: Iterable[float],
        vec1: Iterable[float],
        vec2: Optional[Iterable[float]] = None,
        vec3: Optional[Iterable[float]] = None,
    ) -> None:
        super().__init__(name, "RECTANGULAR", base_point)

        self._vectors: list[Tuple[float, float, float]] = []
        self._vectors.append(self._to_triplet(vec1, "vec1"))

        if vec2 is not None:
            self._vectors.append(self._to_triplet(vec2, "vec2"))
        if vec3 is not None:
            self._vectors.append(self._to_triplet(vec3, "vec3"))

        if len(self._vectors) > 3:
            raise ValueError("Rectangular systems support at most three vectors.")

    @property
    def vectors(self) -> tuple[Tuple[float, float, float], ...]:
        return tuple(self._vectors)

    def to_femaster(self) -> str:
        lines = [self._header(), self._format_triplet(self.base_point)]
        for vec in self._vectors:
            lines.append(self._format_triplet(vec))
        return "\n".join(lines)
