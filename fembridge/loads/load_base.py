from __future__ import annotations

from typing import Callable, Optional, Tuple, Union
from ..coordinates.coordinate_system_base import CoordinateSystem


TIME_DEPENDENT_ERROR = "Time-dependent loads are not supported by FEMaster yet."

AmplitudeFunction = Callable[[float], Union[float, Tuple[float, ...]]]


class Load:
    """Common load behaviour providing optional coordinate system handling."""

    load_type: str = "LOAD"

    def __init__(
        self,
        coordinate_system: Optional[CoordinateSystem] = None,
        amplitude: Optional[AmplitudeFunction] = None,
    ) -> None:
        self.coordinate_system: Optional[CoordinateSystem] = coordinate_system
        self.amplitude: Optional[AmplitudeFunction] = amplitude
        self.is_time_dependent: bool = amplitude is not None

    def _header(self, collector_name: str) -> str:
        header = f"*{self.load_type}, LOAD_COLLECTOR={collector_name}"
        if self.coordinate_system:
            header += f", CSYS={self.coordinate_system.name}"
        return header

    def to_femaster(self, collector_name: str) -> str:  # pragma: no cover - interface only
        raise NotImplementedError("Load subclasses must implement to_femaster().")
