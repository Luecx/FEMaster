"""Repository for reusable load objects and amplitudes."""

from __future__ import annotations

from typing import Iterator

from femaster_api.utils import normalize_name

from .amplitude import Amplitude
from .inertial_load import InertialLoad
from .nodal_force import NodalForce
from .pressure_load import PressureLoad
from .surface_traction import SurfaceTraction
from .thermal_load import ThermalLoad
from .volume_load import VolumeLoad

LoadEntry = NodalForce | SurfaceTraction | PressureLoad | VolumeLoad | ThermalLoad | InertialLoad


class LoadRepository:
    """Repository for reusable load objects and amplitudes."""

    def __init__(self) -> None:
        self._loads: list[LoadEntry] = []
        self._amplitudes: dict[str, Amplitude] = {}

    def add(self, entry: LoadEntry) -> LoadEntry:
        if not isinstance(entry, (NodalForce, SurfaceTraction, PressureLoad, VolumeLoad, ThermalLoad, InertialLoad)):
            raise TypeError("entry must be a load object")
        self._loads.append(entry)
        return entry

    def add_amplitude(self, amplitude: Amplitude) -> Amplitude:
        if not isinstance(amplitude, Amplitude):
            raise TypeError("amplitude must be an Amplitude")
        item = Amplitude(
            normalize_name(amplitude.name),
            tuple((float(time), float(value)) for time, value in amplitude.points),
            amplitude.interpolation,
        )
        self._amplitudes[item.name] = item
        return item

    def amplitudes(self) -> tuple[Amplitude, ...]:
        return tuple(self._amplitudes[key] for key in sorted(self._amplitudes))

    def has_amplitude(self, value: str | Amplitude) -> bool:
        if isinstance(value, Amplitude):
            return any(item is value for item in self._amplitudes.values())
        return normalize_name(value) in self._amplitudes

    def all(self) -> tuple[LoadEntry, ...]:
        return tuple(self._loads)

    def __getitem__(self, index: int | slice) -> LoadEntry | tuple[LoadEntry, ...]:
        if isinstance(index, slice):
            return tuple(self._loads[index])
        return self._loads[index]

    def __iter__(self) -> Iterator[LoadEntry]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._loads)
