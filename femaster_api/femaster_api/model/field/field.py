"""Sparse semantic field."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.named_object import NamedObject
from .field_domain import FieldDomain
from .field_type import FieldType
from .typing import FieldKey


class Field(NamedObject):
    """Sparse field with explicit storage domain, semantics and components."""

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
        if self.components:
            return len(self.components)
        if self.values:
            return len(next(iter(self.values.values())))
        return 0

    def set(self, key: FieldKey, values: Iterable[float]) -> "Field":
        row = tuple(float(value) for value in values)
        if self.components and len(row) != len(self.components):
            raise ValueError(
                f"field {self.name!r} expects {len(self.components)} values, got {len(row)}"
            )
        if self.values and len(row) != self.cols:
            raise ValueError(f"field {self.name!r} has inconsistent row width")
        self.values[key] = row
        return self

    def get(self, key: FieldKey) -> tuple[float, ...]:
        return self.values[key]

    def export(self) -> str:
        lines = [
            keyword(
                "FIELD",
                NAME=self.name,
                TYPE=self.domain.value,
                COLS=self.cols,
                FILL="ZERO",
            )
        ]

        expected_indices = {
            FieldDomain.NODE: 1,
            FieldDomain.ELEMENT: 1,
            FieldDomain.ELEMENT_NODAL: 2,
            FieldDomain.ELEMENT_IP: 2,
            FieldDomain.ELEMENT_MP: 3,
        }.get(self.domain)

        for key, values in self.values.items():
            row_key = key if isinstance(key, tuple) else (key,)
            if expected_indices is not None and len(row_key) != expected_indices:
                raise ValueError(
                    f"field {self.name!r} on {self.domain.value} expects "
                    f"{expected_indices} index columns, got {len(row_key)}"
                )
            lines.append(csv((*row_key, *values)))

        return block(lines)

    def __getitem__(self, key: FieldKey) -> tuple[float, ...]:
        return self.get(key)

    def __len__(self) -> int:
        return len(self.values)
