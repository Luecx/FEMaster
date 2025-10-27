from __future__ import annotations

from typing import Optional, TYPE_CHECKING, Union

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ..coordinates.coordinate_system_base import CoordinateSystem

CoordinateSystemLike = Union["CoordinateSystem", str]


class ConstraintBase:
    """Common base class providing optional coordinate-system handling."""

    constraint_type: str = "CONSTRAINT"

    def __init__(self, coordinate_system: Optional[CoordinateSystemLike] = None) -> None:
        self.coordinate_system: Optional[str] = None
        if coordinate_system is not None:
            self.coordinate_system = self._resolve_coordinate_system(coordinate_system)

    @staticmethod
    def _resolve_coordinate_system(value: CoordinateSystemLike) -> str:
        name = getattr(value, "name", value)  # type: ignore[assignment]
        name_str = str(name).strip()
        if not name_str:
            raise ValueError("Coordinate system reference must have a non-empty name.")
        return name_str

    def to_femaster(self) -> str:  # pragma: no cover - interface only
        raise NotImplementedError("Constraint subclasses must implement to_femaster().")

    def to_asami(self) -> str:  # pragma: no cover - interface only
        raise NotImplementedError("Constraint subclasses must implement to_asami().")
