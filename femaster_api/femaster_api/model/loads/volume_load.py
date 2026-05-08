"""Volume load model object."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from femaster_api.model.elements import Element
from femaster_api.model.loads.amplitude import Amplitude
from femaster_api.model.orientations import Orientation
from femaster_api.model.sets import ElementSet

ElementTarget = Element | ElementSet
Vector3 = tuple[float, float, float]


@dataclass(frozen=True, slots=True, init=False)
class VolumeLoad:
    """Body load applied to an element or element set."""

    target: ElementTarget
    values: Vector3
    orientation: Orientation | None = None
    amplitude: Amplitude | None = None

    def __init__(
        self,
        target: ElementTarget,
        *,
        bx: float | None = None,
        by: float | None = None,
        bz: float | None = None,
        orientation: Orientation | None = None,
        amplitude: Amplitude | None = None,
        values: Iterable[float | None] | None = None,
    ) -> None:
        load_values = _three_values(values if values is not None else (bx, by, bz))
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "values", load_values)
        object.__setattr__(self, "orientation", orientation)
        object.__setattr__(self, "amplitude", amplitude)


def _three_values(values: Iterable[float | None]) -> Vector3:
    result = tuple(0.0 if value is None else float(value) for value in values)
    if len(result) != 3:
        raise ValueError("volume load values must contain three entries")
    return result  # type: ignore[return-value]
