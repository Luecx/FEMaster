"""Sparse field representation for editable model data and solver results.

``Field`` is shared by input-model fields and the post-processing hierarchy, but
the two use cases deliberately have different address semantics.  Editable model
fields store actual ``Node`` / ``Element`` objects.  Only element-local indices
remain integers because they are intrinsic positions inside an element rather
than references to another model object.

RES/FRD result fields do not own the originating model graph and therefore keep
the semantic solver addresses encoded by the result file.  Result readers use a
private insertion path for those serialized addresses; the public ``set`` method
is reserved for object-valued editable-model addressing.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject
from ..element.element import Element
from ..node.node import Node
from .field_domain import FieldDomain
from .field_type import FieldType
from .typing import FieldKey, ModelFieldKey, ResultFieldKey


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

    def set(self, key: ModelFieldKey, values: Iterable[float]) -> "Field":
        """Store one editable-model row addressed by concrete model objects.

        ``Node`` and ``Element`` IDs are intentionally not accepted here.  For
        element-local domains the first tuple entry is the concrete ``Element``;
        subsequent entries are zero-based local node/IP/material-point indices.
        """

        self._validate_model_key(key)
        self._store(key, values)
        return self

    def _set_result(
        self,
        key: ResultFieldKey,
        values: Iterable[float],
    ) -> "Field":
        """Store one serialized result row without pretending it is a model link."""

        self._store(key, values)
        return self

    def _store(self, key: FieldKey, values: Iterable[float]) -> None:
        """Store one row after validating component width."""

        row = tuple(float(value) for value in values)

        if self.components and len(row) != len(self.components):
            raise ValueError(
                f"field {self.name!r} expects {len(self.components)} values, "
                f"got {len(row)}"
            )

        if self.values and len(row) != self.cols:
            raise ValueError(f"field {self.name!r} has inconsistent row width")

        self.values[key] = row

    def get(self, key: FieldKey) -> tuple[float, ...]:
        """Return one field row by its model object or result address."""

        return self.values[key]

    def export(self) -> str:
        """Export an editable model field as one native ``*FIELD`` block."""

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
            address = self._model_address(key)
            lines.append(csv((*address, *self.values[key])))

        return block(lines)

    def _validate_model_key(self, key: ModelFieldKey) -> None:
        """Validate object-valued addressing for the selected model domain."""

        if self.domain is FieldDomain.NODE:
            if not isinstance(key, Node):
                raise TypeError("NODE field keys must be Node objects")
            return

        if self.domain is FieldDomain.ELEMENT:
            if not isinstance(key, Element):
                raise TypeError("ELEMENT field keys must be Element objects")
            return

        if self.domain in {FieldDomain.ELEMENT_NODAL, FieldDomain.ELEMENT_IP}:
            if not (
                isinstance(key, tuple)
                and len(key) == 2
                and isinstance(key[0], Element)
                and type(key[1]) is int
            ):
                raise TypeError(
                    f"{self.domain.value} field keys must be (Element, int)"
                )
            if key[1] < 0:
                raise ValueError("element-local field indices must be non-negative")
            return

        if self.domain is FieldDomain.ELEMENT_MP:
            if not (
                isinstance(key, tuple)
                and len(key) == 3
                and isinstance(key[0], Element)
                and type(key[1]) is int
                and type(key[2]) is int
            ):
                raise TypeError(
                    "ELEMENT_MP field keys must be (Element, int, int)"
                )
            if key[1] < 0 or key[2] < 0:
                raise ValueError("element-local field indices must be non-negative")
            return

        raise ValueError(
            f"field {self.name!r} with domain {self.domain.value} "
            "has no editable-model address representation"
        )

    def _model_address(self, key: FieldKey) -> tuple[int, ...]:
        """Reduce a validated model-object key to native numeric address columns."""

        self._validate_model_key(key)  # type: ignore[arg-type]

        if isinstance(key, Node):
            return (key.id,)
        if isinstance(key, Element):
            return (key.id,)

        element = key[0]
        return (element.id, *key[1:])  # type: ignore[union-attr]

    @staticmethod
    def _sort_key(key: FieldKey) -> tuple[str, ...]:
        """Build a deterministic ordering for model and result field addresses."""

        values = key if isinstance(key, tuple) else (key,)
        normalized: list[str] = []
        for value in values:
            if isinstance(value, (Node, Element)):
                normalized.append(str(value.id))
            else:
                normalized.append(str(value))
        return tuple(normalized)

    def __getitem__(self, key: FieldKey) -> tuple[float, ...]:
        return self.get(key)

    def __len__(self) -> int:
        return len(self.values)
