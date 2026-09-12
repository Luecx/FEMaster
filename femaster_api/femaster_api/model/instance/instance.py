"""Rigid assembly instance of one reusable ``Part`` object.

An ``Instance`` stores the concrete ``Part`` it instantiates; a part name is not
accepted as a substitute object.  The native ``PART=...`` token is derived from
``part.name`` only during export.  This guarantees that an instance cannot carry
a dangling part reference inside an otherwise valid Python project.

Optional translation and axis-angle rotation rows mirror FEMaster's native
``*INSTANCE`` syntax and are validated at construction time.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject
from ..part.part import Part


class Instance(NamedObject):
    """Rigid placement of one concrete ``Part`` inside the assembly."""

    def __init__(
        self,
        name: str,
        part: Part,
        *,
        translation: Iterable[float] | None = None,
        rotation: tuple[Iterable[float], Iterable[float], float] | None = None,
    ) -> None:
        super().__init__(name)
        if not isinstance(part, Part):
            raise TypeError("part must be a Part object")
        self.part = part

        self.translation: tuple[float, float, float] | None = None
        if translation is not None:
            values = tuple(float(value) for value in translation)
            if len(values) != 3:
                raise ValueError("instance translation requires exactly 3 values")
            self.translation = values

        self.rotation: tuple[
            tuple[float, float, float],
            tuple[float, float, float],
            float,
        ] | None = None
        if rotation is not None:
            point_a = tuple(float(value) for value in rotation[0])
            point_b = tuple(float(value) for value in rotation[1])
            if len(point_a) != 3 or len(point_b) != 3:
                raise ValueError(
                    "instance rotation axis points require exactly 3 values"
                )
            self.rotation = (point_a, point_b, float(rotation[2]))

    def export(self) -> str:
        """Export one native ``*INSTANCE`` block including optional placement."""

        lines = [keyword("INSTANCE", NAME=self.name, PART=self.part.name)]
        if self.translation is not None:
            lines.append(csv(self.translation))
        if self.rotation is not None:
            point_a, point_b, angle = self.rotation
            lines.append(csv((*point_a, *point_b, angle)))
        return block(lines)
