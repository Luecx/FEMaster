from __future__ import annotations

from typing import Iterable, Tuple


class CoordinateSystem:
    """Base class for all coordinate system representations."""

    def __init__(self, name: str, system_type: str, base_point: Iterable[float]):
        if not name:
            raise ValueError("Coordinate system name must be non-empty.")
        if not system_type:
            raise ValueError("Coordinate system type must be non-empty.")

        self.name: str = str(name)
        self.system_type: str = str(system_type)
        self.base_point: Tuple[float, float, float] = self._to_triplet(base_point, "base_point")

    @staticmethod
    def _to_triplet(values: Iterable[float], label: str) -> Tuple[float, float, float]:
        try:
            triple = tuple(float(v) for v in values)
        except TypeError as exc:
            raise TypeError(f"{label} must be an iterable of three numeric values.") from exc
        if len(triple) != 3:
            raise ValueError(f"{label} must contain exactly three numeric values.")
        return triple

    @staticmethod
    def _format_triplet(values: Tuple[float, float, float]) -> str:
        return f"{values[0]}, {values[1]}, {values[2]}"

    def _header(self) -> str:
        return f"*COORDINATE_SYSTEM, NAME={self.name}, TYPE={self.system_type}"

    def to_femaster(self) -> str:
        raise NotImplementedError("CoordinateSystem subclasses must implement to_femaster().")

    def to_asami(self) -> str:
        raise NotImplementedError("CoordinateSystem.to_asami not implemented yet.")
