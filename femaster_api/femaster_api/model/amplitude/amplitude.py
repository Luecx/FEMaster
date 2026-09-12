"""Named scalar amplitude."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject
from .amplitude_interpolation import AmplitudeInterpolation


class Amplitude(NamedObject):
    """Named scalar time function shared by loads and future BCs."""

    def __init__(
        self,
        name: str,
        points: Iterable[tuple[float, float]] = (),
        interpolation: AmplitudeInterpolation = AmplitudeInterpolation.LINEAR,
    ) -> None:
        super().__init__(name)
        self.points = [(float(time), float(value)) for time, value in points]
        self.interpolation = interpolation

    def add(self, time: float, value: float) -> "Amplitude":
        self.points.append((float(time), float(value)))
        return self

    def export(self) -> str:
        return block([
            keyword("AMPLITUDE", NAME=self.name, TYPE=self.interpolation.value),
            *(csv(point) for point in self.points),
        ])
