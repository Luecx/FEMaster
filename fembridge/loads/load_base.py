from __future__ import annotations

from typing import Optional, Union
from ..coordinates.coordinate_system_base import CoordinateSystem



class Load:
    """Common load behaviour providing optional coordinate system handling."""

    load_type: str = "LOAD"

    def __init__(self, coordinate_system: Optional[CoordinateSystem] = None) -> None:
        self.coordinate_system: Optional[str] = None
        if coordinate_system is not None:
            self.coordinate_system = self._resolve_coordinate_system(coordinate_system)

    @staticmethod
    def _resolve_coordinate_system(value: CoordinateSystem) -> str:
        if hasattr(value, "name"):
            name = getattr(value, "name")
        else:
            name = value  # type: ignore[assignment]
        name_str = str(name).strip()
        if not name_str:
            raise ValueError("Coordinate system reference must have a non-empty name.")
        return name_str

    def _header(self, collector_name: str) -> str:
        header = f"*{self.load_type}, LOAD_COLLECTOR={collector_name}"
        if self.coordinate_system:
            header += f", CSYS={self.coordinate_system}"
        return header

    def to_femaster(self, collector_name: str) -> str:  # pragma: no cover - interface only
        raise NotImplementedError("Load subclasses must implement to_femaster().")
