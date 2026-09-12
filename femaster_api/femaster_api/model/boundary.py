"""FEMaster amplitudes, loads and prescribed supports.

Loads and supports are owned directly by named collectors. There is no second
global repository containing the same entries. Analysis steps activate collector
names, matching FEMaster's load/support collector semantics.
"""

from __future__ import annotations

from enum import Enum
from typing import Iterable

from .._format import block, csv, keyword
from ..repository import NamedObject, NamedRepository
from .typing import EntityReference


class AmplitudeInterpolation(Enum):
    """Interpolation mode between amplitude samples."""

    LINEAR = "LINEAR"
    STEP   = "STEP"


class Amplitude(NamedObject):
    """Named scalar time function shared by load definitions."""

    def __init__(
        self,
        name: str,
        points: Iterable[tuple[float, float]] = (),
        interpolation: AmplitudeInterpolation = AmplitudeInterpolation.LINEAR,
    ) -> None:
        super().__init__(name)
        self.points        = [(float(time), float(value)) for time, value in points]
        self.interpolation = interpolation

    def add(self, time: float, value: float) -> "Amplitude":
        """Append one sample and return this amplitude."""

        self.points.append((float(time), float(value)))
        return self

    def to_femaster(self) -> str:
        return block([
            keyword("AMPLITUDE", NAME=self.name, TYPE=self.interpolation.value),
            *(csv(point) for point in self.points),
        ])


class AmplitudeRepository(NamedRepository[Amplitude]):
    """Repository of named global amplitudes."""

    def to_femaster(self) -> str:
        return "\n\n".join(amplitude.to_femaster() for amplitude in self)


class Load:
    """Base class for one collector-owned FEMaster load definition."""

    def to_femaster(self, collector: str) -> str:
        raise NotImplementedError


class NodalForce(Load):
    """Concentrated nodal force and moment applied to a node or node region."""

    def __init__(
        self,
        target: EntityReference,
        values: Iterable[float] = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        *,
        orientation: str | None = None,
        amplitude: str | None = None,
    ) -> None:
        self.target      = target
        self.values      = tuple(float(value) for value in values)
        self.orientation = orientation
        self.amplitude   = amplitude
        if len(self.values) != 6:
            raise ValueError("NodalForce requires exactly 6 force/moment values")

    def to_femaster(self, collector: str) -> str:
        return block([
            keyword(
                "CLOAD",
                LOAD_COLLECTOR=collector,
                ORIENTATION=self.orientation,
                AMPLITUDE=self.amplitude,
            ),
            csv((self.target, *self.values)),
        ])


class SurfaceTraction(Load):
    """Distributed traction vector applied to a surface region."""

    def __init__(
        self,
        target: EntityReference,
        values: Iterable[float],
        *,
        orientation: str | None = None,
        amplitude: str | None = None,
    ) -> None:
        self.target      = target
        self.values      = tuple(float(value) for value in values)
        self.orientation = orientation
        self.amplitude   = amplitude
        if len(self.values) != 3:
            raise ValueError("SurfaceTraction requires exactly 3 components")

    def to_femaster(self, collector: str) -> str:
        return block([
            keyword(
                "DLOAD",
                LOAD_COLLECTOR=collector,
                ORIENTATION=self.orientation,
                AMPLITUDE=self.amplitude,
            ),
            csv((self.target, *self.values)),
        ])


class PressureLoad(Load):
    """Scalar pressure applied to a surface region."""

    def __init__(self, target: EntityReference, pressure: float, *, amplitude: str | None = None) -> None:
        self.target    = target
        self.pressure  = float(pressure)
        self.amplitude = amplitude

    def to_femaster(self, collector: str) -> str:
        return block([
            keyword("PLOAD", LOAD_COLLECTOR=collector, AMPLITUDE=self.amplitude),
            csv((self.target, self.pressure)),
        ])


class VolumeLoad(Load):
    """Distributed body-force vector applied to an element region."""

    def __init__(
        self,
        target: EntityReference,
        values: Iterable[float],
        *,
        orientation: str | None = None,
        amplitude: str | None = None,
    ) -> None:
        self.target      = target
        self.values      = tuple(float(value) for value in values)
        self.orientation = orientation
        self.amplitude   = amplitude
        if len(self.values) != 3:
            raise ValueError("VolumeLoad requires exactly 3 components")

    def to_femaster(self, collector: str) -> str:
        return block([
            keyword(
                "VLOAD",
                LOAD_COLLECTOR=collector,
                ORIENTATION=self.orientation,
                AMPLITUDE=self.amplitude,
            ),
            csv((self.target, *self.values)),
        ])


class ThermalLoad(Load):
    """Structural thermal load referencing a named nodal temperature field."""

    def __init__(self, temperature_field: str, reference_temperature: float = 0.0) -> None:
        self.temperature_field      = str(temperature_field)
        self.reference_temperature = float(reference_temperature)

    def to_femaster(self, collector: str) -> str:
        return keyword(
            "TLOAD",
            LOAD_COLLECTOR=collector,
            TEMPERATUREFIELD=self.temperature_field,
            REFERENCETEMPERATURE=self.reference_temperature,
        )


class InertialLoad(Load):
    """Rigid-body translational and rotational inertia loading."""

    def __init__(
        self,
        target: str,
        *,
        center: Iterable[float] = (0.0, 0.0, 0.0),
        center_acceleration: Iterable[float] = (0.0, 0.0, 0.0),
        omega: Iterable[float] = (0.0, 0.0, 0.0),
        alpha: Iterable[float] = (0.0, 0.0, 0.0),
        consider_point_masses: bool = True,
    ) -> None:
        self.target                = str(target)
        self.center                = tuple(float(value) for value in center)
        self.center_acceleration   = tuple(float(value) for value in center_acceleration)
        self.omega                 = tuple(float(value) for value in omega)
        self.alpha                 = tuple(float(value) for value in alpha)
        self.consider_point_masses = bool(consider_point_masses)

    def to_femaster(self, collector: str) -> str:
        return block([
            keyword(
                "INERTIALOAD",
                LOAD_COLLECTOR=collector,
                CONSIDER_POINT_MASSES=int(self.consider_point_masses),
            ),
            csv((
                self.target,
                *self.center,
                *self.center_acceleration,
                *self.omega,
                *self.alpha,
            )),
        ])


class LoadCollector(NamedObject):
    """Named collection of loads activated together by analysis steps."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.loads: list[Load] = []

    def add(self, load: Load) -> Load:
        if not isinstance(load, Load):
            raise TypeError("load must derive from Load")
        self.loads.append(load)
        return load

    def to_femaster(self) -> str:
        return "\n\n".join(load.to_femaster(self.name) for load in self.loads)


class LoadCollectorRepository(NamedRepository[LoadCollector]):
    """Repository of named load collectors."""

    def to_femaster(self) -> str:
        return "\n\n".join(collector.to_femaster() for collector in self if collector.loads)


class Support:
    """Prescribed structural DOFs on a node or node region.

    ``None`` leaves a DOF unconstrained. Numeric values prescribe the matching
    translational or rotational DOF in FEMaster's six-component order.
    """

    def __init__(
        self,
        target: EntityReference,
        values: Iterable[float | None],
        *,
        orientation: str | None = None,
    ) -> None:
        self.target      = target
        self.values      = tuple(None if value is None else float(value) for value in values)
        self.orientation = orientation
        if len(self.values) > 6:
            raise ValueError("Support accepts at most 6 prescribed DOF values")

    def to_femaster(self, collector: str) -> str:
        values = list(self.values)
        while values and values[-1] is None:
            values.pop()
        return block([
            keyword(
                "SUPPORT",
                SUPPORT_COLLECTOR=collector,
                ORIENTATION=self.orientation,
            ),
            csv((self.target, *values)),
        ])


class SupportCollector(NamedObject):
    """Named collection of prescribed supports activated by analysis steps."""

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.supports: list[Support] = []

    def add(self, support: Support) -> Support:
        if not isinstance(support, Support):
            raise TypeError("support must be a Support")
        self.supports.append(support)
        return support

    def to_femaster(self) -> str:
        return "\n\n".join(support.to_femaster(self.name) for support in self.supports)


class SupportCollectorRepository(NamedRepository[SupportCollector]):
    """Repository of named support collectors."""

    def to_femaster(self) -> str:
        return "\n\n".join(collector.to_femaster() for collector in self if collector.supports)
