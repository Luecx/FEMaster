"""Object-keyed input field data object."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Iterable

from .field_domain import FieldDomain


@dataclass(frozen=True, slots=True)
class Field:
    """Named tabular input field keyed by model objects.

    Node fields should use ``Node`` keys, element fields should use
    ``Element`` keys, and integration-point fields may use explicit tuple keys
    such as ``(element, integration_point)``.
    """

    name: str
    domain: FieldDomain
    cols: int
    values: dict[object, tuple[float | None, ...]] = field(default_factory=dict)
    fill: str = "ZERO"

    def set(self, key: object, values: Iterable[float | None]) -> "Field":
        """Return a copy with one object-keyed field row set."""

        if key is None:
            raise ValueError("field key cannot be None")
        row_values = tuple(None if value is None else float(value) for value in values)
        if len(row_values) > self.cols:
            raise ValueError(f"field {self.name!r} has {self.cols} columns, got {len(row_values)} values")
        padded = row_values + (None,) * (self.cols - len(row_values))
        new_values = dict(self.values)
        new_values[key] = padded
        return replace(self, values=new_values)


InputField = Field
