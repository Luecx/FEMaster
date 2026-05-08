"""Surface traction model object."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from femaster_api.model.loads.amplitude import Amplitude
from femaster_api.model.orientations import Orientation
from femaster_api.model.sets import SurfaceSet
from femaster_api.model.surfaces import SurfaceDefinition

SurfaceTarget = SurfaceDefinition | SurfaceSet
Vector3 = tuple[float, float, float]


@dataclass(frozen=True, slots=True, init=False)
class SurfaceTraction:
    """Distributed traction applied to a surface or surface set."""

    target: SurfaceTarget
    values: Vector3
    orientation: Orientation | None = None
    amplitude: Amplitude | None = None

    def __init__(
        self,
        target: SurfaceTarget,
        *,
        tx: float | None = None,
        ty: float | None = None,
        tz: float | None = None,
        orientation: Orientation | None = None,
        amplitude: Amplitude | None = None,
        values: Iterable[float | None] | None = None,
    ) -> None:
        traction_values = _three_values(values if values is not None else (tx, ty, tz))
        object.__setattr__(self, "target", target)
        object.__setattr__(self, "values", traction_values)
        object.__setattr__(self, "orientation", orientation)
        object.__setattr__(self, "amplitude", amplitude)


def _three_values(values: Iterable[float | None]) -> Vector3:
    result = tuple(0.0 if value is None else float(value) for value in values)
    if len(result) != 3:
        raise ValueError("surface traction values must contain three entries")
    return result  # type: ignore[return-value]
