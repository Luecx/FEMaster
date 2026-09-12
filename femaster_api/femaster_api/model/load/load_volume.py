"""Distributed body-force vector on concrete element-domain objects.

``VolumeLoad`` accepts an ``Element`` or ``ElementRegion`` target.  Orientation
and amplitude are likewise concrete shared model objects.  IDs and semantic names
are emitted only by ``export`` so the Python graph never stores those tokens as
surrogates for relationships.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..amplitude.amplitude import Amplitude
from ..common.format import block, csv, keyword
from ..coordinate_system.coordinate_system import CoordinateSystem
from ..element.element import Element
from ..region.region_element import ElementRegion
from .load import Load


class VolumeLoad(Load):
    """Three-component body force applied to an element-domain target."""

    def __init__(
        self,
        target: Element | ElementRegion,
        values: Iterable[float],
        *,
        orientation: CoordinateSystem | None = None,
        amplitude: Amplitude | None = None,
    ) -> None:
        if not isinstance(target, (Element, ElementRegion)):
            raise TypeError("target must be Element or ElementRegion")
        if orientation is not None and not isinstance(orientation, CoordinateSystem):
            raise TypeError("orientation must be a CoordinateSystem object")
        if amplitude is not None and not isinstance(amplitude, Amplitude):
            raise TypeError("amplitude must be an Amplitude object")

        self.target = target
        self.values = tuple(float(value) for value in values)
        self.orientation = orientation
        self.amplitude = amplitude

        if len(self.values) != 3:
            raise ValueError("VolumeLoad requires exactly 3 components")

    def export(self, collector: str) -> str:
        target = self.target.id if isinstance(self.target, Element) else self.target.name
        return block([
            keyword(
                "VLOAD",
                LOAD_COLLECTOR=collector,
                ORIENTATION=(
                    self.orientation.name if self.orientation is not None else None
                ),
                AMPLITUDE=self.amplitude.name if self.amplitude is not None else None,
            ),
            csv((target, *self.values)),
        ])
