"""Rigid-body translational and rotational inertia loading on element objects.

The load stores an ``Element`` or ``ElementRegion`` target together with center
position, center acceleration, angular velocity and angular acceleration vectors.
The target's ID/name is a file-format token only and is never retained as the
in-memory relationship.

All vector dimensions are validated at construction time so export cannot
silently emit malformed ``*INERTIALOAD`` rows.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..element.element import Element
from ..region.region_element import ElementRegion
from .load import Load


class InertialLoad(Load):
    """Combined translational and rotational inertia loading."""

    def __init__(
        self,
        target: Element | ElementRegion,
        *,
        center: Iterable[float] = (0.0, 0.0, 0.0),
        center_acceleration: Iterable[float] = (0.0, 0.0, 0.0),
        omega: Iterable[float] = (0.0, 0.0, 0.0),
        alpha: Iterable[float] = (0.0, 0.0, 0.0),
        consider_point_masses: bool = True,
    ) -> None:
        if not isinstance(target, (Element, ElementRegion)):
            raise TypeError("target must be Element or ElementRegion")
        self.target = target
        self.center = self._vec3(center, "center")
        self.center_acceleration = self._vec3(
            center_acceleration,
            "center_acceleration",
        )
        self.omega = self._vec3(omega, "omega")
        self.alpha = self._vec3(alpha, "alpha")
        self.consider_point_masses = bool(consider_point_masses)

    @staticmethod
    def _vec3(
        values: Iterable[float],
        name: str,
    ) -> tuple[float, float, float]:
        result = tuple(float(value) for value in values)
        if len(result) != 3:
            raise ValueError(f"{name} requires exactly 3 components")
        return result

    def export(self, collector: str) -> str:
        target = self.target.id if isinstance(self.target, Element) else self.target.name
        return block([
            keyword(
                "INERTIALOAD",
                LOAD_COLLECTOR=collector,
                CONSIDER_POINT_MASSES=int(self.consider_point_masses),
            ),
            csv((
                target,
                *self.center,
                *self.center_acceleration,
                *self.omega,
                *self.alpha,
            )),
        ])
