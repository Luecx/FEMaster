"""Named reusable collection of loads."""

from __future__ import annotations

from dataclasses import dataclass, replace

from .inertial_load import InertialLoad
from .nodal_force import NodalForce
from .pressure_load import PressureLoad
from .surface_traction import SurfaceTraction
from .thermal_load import ThermalLoad
from .volume_load import VolumeLoad

LoadEntry = NodalForce | SurfaceTraction | PressureLoad | VolumeLoad | ThermalLoad | InertialLoad


@dataclass(frozen=True, slots=True)
class LoadCollector:
    """Named list of loads that can be attached to loadcases."""

    name: str
    loads: tuple[LoadEntry, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "loads", tuple(self.loads))

    def add(self, load: LoadEntry) -> "LoadCollector":
        return replace(self, loads=(*self.loads, load))
