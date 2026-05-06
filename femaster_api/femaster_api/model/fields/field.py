"""Generic tabular field data object."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Iterable

from .field_domain import FieldDomain


@dataclass(frozen=True, slots=True)
class Field:
    name: str
    domain: FieldDomain
    cols: int
    values: dict[int, tuple[float | None, ...]] = field(default_factory=dict)
    fill: str = "ZERO"

    def set(self, id: int, values: Iterable[float | None]) -> "Field":
        if not isinstance(id, int) or id < 0:
            raise ValueError("field row id must be a non-negative integer")
        row_id = id
        row_values = tuple(None if value is None else float(value) for value in values)
        if len(row_values) > self.cols:
            raise ValueError(f"field {self.name!r} has {self.cols} columns, got {len(row_values)} values")
        padded = row_values + (None,) * (self.cols - len(row_values))
        new_values = dict(self.values)
        new_values[row_id] = padded
        return replace(self, values=new_values)
