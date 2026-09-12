"""Named scalar amplitude shared by time-dependent model definitions.

An amplitude stores ordered ``(time, value)`` samples and the interpolation rule
used between them.  Loads retain the ``Amplitude`` object itself when they use a
time function; ``amplitude.name`` is consulted only while native keywords are
serialized.  The same curve can therefore be shared by identity without copying
data or passing semantic-name strings through the Python model.

The class owns the native ``*AMPLITUDE`` representation directly.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject
from .amplitude_interpolation import AmplitudeInterpolation


class Amplitude(NamedObject):
    """Named scalar time function exported through ``*AMPLITUDE``."""

    def __init__(
        self,
        name: str,
        points: Iterable[tuple[float, float]] = (),
        interpolation: AmplitudeInterpolation = AmplitudeInterpolation.LINEAR,
    ) -> None:
        super().__init__(name)
        self.points = [
            (float(time), float(value))
            for time, value in points
        ]
        self.interpolation = interpolation

    def add(self, time: float, value: float) -> "Amplitude":
        """Append one sample and return this amplitude."""

        self.points.append((float(time), float(value)))
        return self

    def export(self) -> str:
        """Export one native ``*AMPLITUDE`` block."""

        return block([
            keyword(
                "AMPLITUDE",
                NAME=self.name,
                TYPE=self.interpolation.value,
            ),
            *(csv(point) for point in self.points),
        ])
