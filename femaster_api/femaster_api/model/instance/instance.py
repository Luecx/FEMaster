"""Rigid Part instance placement."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject


class Instance(NamedObject):
    """Rigid placement of one reusable Part in the assembly."""

    def __init__(
        self,
        name: str,
        part: str,
        *,
        translation: Iterable[float] | None = None,
        rotation: tuple[Iterable[float], Iterable[float], float] | None = None,
    ) -> None:
        super().__init__(name)
        self.part = str(part)

        self.translation = None
        if translation is not None:
            values = tuple(float(value) for value in translation)
            if len(values) != 3:
                raise ValueError("instance translation requires exactly 3 values")
            self.translation = values

        self.rotation = None
        if rotation is not None:
            point_a = tuple(float(value) for value in rotation[0])
            point_b = tuple(float(value) for value in rotation[1])
            if len(point_a) != 3 or len(point_b) != 3:
                raise ValueError("instance rotation axis points require exactly 3 values")
            self.rotation = (point_a, point_b, float(rotation[2]))

    def export(self) -> str:
        lines = [keyword("INSTANCE", NAME=self.name, PART=self.part)]
        if self.translation is not None:
            lines.append(csv(self.translation))
        if self.rotation is not None:
            point_a, point_b, angle = self.rotation
            lines.append(csv((*point_a, *point_b, angle)))
        return block(lines)
