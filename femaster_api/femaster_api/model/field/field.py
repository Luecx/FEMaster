"""Sparse semantic field shared by model input and solver results.

``Field`` combines three independent pieces of information: a human/native name,
a physical ``FieldDomain`` and a canonical ``FieldType``.  Values are stored by
semantic entity identifiers rather than dense solver row numbers.  This allows
native RES identifiers such as ``17`` and ``"bolt.17"`` to survive a read/write
cycle without global renumbering.

For model input export the addressing shape is validated against the domain:
NODE/ELEMENT use one identifier, ELEMENT_NODAL/ELEMENT_IP use an element plus
one local index, and ELEMENT_MP uses element, local integration point and local
material point.  UNKNOWN result fields are intentionally not exportable as
``*FIELD`` because the input parser has no UNKNOWN domain.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject
from .field_domain import FieldDomain
from .field_type import FieldType
from .typing import FieldKey


class Field(NamedObject):
    """Sparse field with explicit domain, semantics and component names."""

    def __init__(
        self,
        name: str,
        domain: FieldDomain,
        components: Iterable[str] = (),
        *,
        type: FieldType | None = None,
    ) -> None:
        super().__init__(name)
        self.domain = domain
        self.type = type or FieldType.from_name(name)
        self.components = tuple(str(component) for component in components)
        self.values: dict[FieldKey, tuple[float, ...]] = {}

    @property
    def cols(self) -> int:
        """Return the numerical component count of one row."""

        if self.components:
            return len(self.components)
        if self.values:
            return len(next(iter(self.values.values())))
        return 0

    def set(self, key: FieldKey, values: Iterable[float]) -> "Field":
        """Store one semantic row while enforcing a consistent row width."""

        row = tuple(float(value) for value in values)

        if self.components and len(row) != len(self.components):
            raise ValueError(
                f"field {self.name!r} expects {len(self.components)} values, "
                f"got {len(row)}"
            )

        if self.values and len(row) != self.cols:
            raise ValueError(f"field {self.name!r} has inconsistent row width")

        self.values[key] = row
        return self

    def get(self, key: FieldKey) -> tuple[float, ...]:
        """Return one field row by semantic key."""

        return self.values[key]

    def export(self) -> str:
        """Export this field as a native FEMaster ``*FIELD`` block."""

        if self.domain is FieldDomain.UNKNOWN:
            raise ValueError(
                f"field {self.name!r} has UNKNOWN domain and cannot be exported"
            )

        lines = [
            keyword(
                "FIELD",
                NAME=self.name,
                TYPE=self.domain.value,
                COLS=self.cols,
                FILL="ZERO",
            )
        ]

        for key in sorted(self.values, key=self._sort_key):
            address = key if isinstance(key, tuple) else (key,)
            self._validate_address(address)
            lines.append(csv((*address, *self.values[key])))

        return block(lines)

    def _validate_address(self, address: tuple[object, ...]) -> None:
        """Validate the number of semantic index columns for this domain."""

        required = {
            FieldDomain.NODE: 1,
            FieldDomain.ELEMENT: 1,
            FieldDomain.ELEMENT_NODAL: 2,
            FieldDomain.ELEMENT_IP: 2,
            FieldDomain.ELEMENT_MP: 3,
        }[self.domain]

        if len(address) != required:
            raise ValueError(
                f"field {self.name!r} with domain {self.domain.value} "
                f"requires {required} index columns, got {len(address)}"
            )

    @staticmethod
    def _sort_key(key: FieldKey) -> tuple[str, ...]:
        values = key if isinstance(key, tuple) else (key,)
        return tuple(str(value) for value in values)

    def __getitem__(self, key: FieldKey) -> tuple[float, ...]:
        return self.get(key)

    def __len__(self) -> int:
        return len(self.values)
