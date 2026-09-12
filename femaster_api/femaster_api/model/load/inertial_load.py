"""Rigid-body inertial loading."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .load import Load


class InertialLoad(Load):
    """Translational and rotational rigid-body inertia loading."""

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
        self.target = str(target)
        self.center = tuple(float(value) for value in center)
        self.center_acceleration = tuple(float(value) for value in center_acceleration)
        self.omega = tuple(float(value) for value in omega)
        self.alpha = tuple(float(value) for value in alpha)
        self.consider_point_masses = bool(consider_point_masses)

        for name, values in (
            ("center", self.center),
            ("center_acceleration", self.center_acceleration),
            ("omega", self.omega),
            ("alpha", self.alpha),
        ):
            if len(values) != 3:
                raise ValueError(f"{name} requires exactly 3 values")

    def export(self, collector: str) -> str:
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
