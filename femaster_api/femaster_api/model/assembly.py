"""Reusable parts and rigid assembly instances.

The default part is not stored separately on Project. It is the permanent object
at ``project.parts[0]`` and ``project.parts.default()`` is exactly that lookup.
Deleting or replacing the default part is forbidden. Explicit parts begin at
repository position one.
"""

from __future__ import annotations

from typing import Iterable, Iterator

from .._format import block, blocks, csv, keyword
from ..repository import NamedObject, NamedRepository
from .mesh import ElementRepository, NodeRepository, Surface
from .regions import RegionRepository
from .sections import SectionRepository


class Part(NamedObject):
    """Reusable part-local finite-element definition.

    A Part owns local mesh topology, named regions, boundary surfaces and section
    assignments. Materials, loads, supports, constraints and steps remain global
    Project definitions.
    """

    def __init__(self, name: str) -> None:
        super().__init__(name)
        self.nodes    = NodeRepository()
        self.elements = ElementRepository()
        self.regions  = RegionRepository()
        self.surfaces = NamedRepository[Surface]()
        self.sections = SectionRepository()

    def to_femaster(self, *, root: bool = False) -> str:
        """Return this part in root scope or as an explicit PART block."""

        body = blocks((
            self.nodes.to_femaster(),
            self.elements.to_femaster(),
            self.regions.to_femaster(),
            "\n\n".join(surface.to_femaster() for surface in self.surfaces),
            self.sections.to_femaster(),
        ))

        # The implicit default part is represented by ordinary root definitions
        # and therefore has no PART/ENDPART wrapper in the input deck.
        if root:
            return body

        lines = [keyword("PART", NAME=self.name)]
        if body:
            lines.append(body)
        lines.append(keyword("ENDPART"))
        return block(lines)


class PartRepository(NamedRepository[Part]):
    """Repository of reusable Parts with a permanent default at position zero."""

    DEFAULT_NAME = "__DEFAULT_PART__"

    def __init__(self) -> None:
        super().__init__()
        super().add(Part(self.DEFAULT_NAME))

    def default(self) -> Part:
        """Return the implicit default part at repository position zero."""

        return self[0]

    def explicit(self) -> Iterator[Part]:
        """Iterate explicit user Parts, excluding the permanent default part."""

        return iter(self._items[1:])

    def remove(self, key: int | str) -> Part:
        """Remove an explicit Part while protecting the default Part."""

        part = self[key]
        if part is self.default():
            raise ValueError("the default part cannot be removed")
        return super().remove(key)

    def clear(self) -> None:
        """Remove all explicit Parts while preserving the default Part."""

        default = self.default()
        self._items = [default]
        self._names = {default.name: default}

    def to_femaster(self) -> str:
        """Return root topology followed by all explicit PART blocks."""

        return blocks((
            self.default().to_femaster(root=True),
            *(part.to_femaster() for part in self.explicit()),
        ))


class Instance(NamedObject):
    """Rigid placement of one reusable Part in the assembled model."""

    def __init__(
        self,
        name: str,
        part: str,
        *,
        translation: Iterable[float] | None = None,
        rotation: tuple[
            Iterable[float],
            Iterable[float],
            float,
        ] | None = None,
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

    def to_femaster(self) -> str:
        """Return the FEMaster INSTANCE block including optional placement rows."""

        lines = [keyword("INSTANCE", NAME=self.name, PART=self.part)]

        if self.translation is not None:
            lines.append(csv(self.translation))

        if self.rotation is not None:
            point_a, point_b, angle = self.rotation
            lines.append(csv((*point_a, *point_b, angle)))

        return block(lines)


class InstanceRepository(NamedRepository[Instance]):
    """Repository of named assembly instances.

    The repository exports only INSTANCE definitions. Project owns the ASSEMBLY
    scope because assembly-level sets and surfaces must share that same scope.
    """

    def to_femaster(self) -> str:
        return "\n\n".join(instance.to_femaster() for instance in self)
