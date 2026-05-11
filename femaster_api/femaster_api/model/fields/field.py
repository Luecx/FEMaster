"""Generic tabular field data object."""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Iterable

from .field_domain import FieldDomain

FieldKey = int | tuple[int, ...]


@dataclass(frozen=True, slots=True)
class Field:
    name: str
    domain: FieldDomain
    cols: int
    values: dict[FieldKey, tuple[float | None, ...]] = field(default_factory=dict)
    fill: str = "ZERO"

    def set(self, id: FieldKey, values: Iterable[float | None]) -> "Field":
        row_id = _normalize_key(id)
        row_values = tuple(None if value is None else float(value) for value in values)
        if len(row_values) > self.cols:
            raise ValueError(f"field {self.name!r} has {self.cols} columns, got {len(row_values)} values")
        padded = row_values + (None,) * (self.cols - len(row_values))
        new_values = dict(self.values)
        new_values[row_id] = padded
        return replace(self, values=new_values)

    def row(self, id: FieldKey) -> tuple[float | None, ...]:
        return self.values[_normalize_key(id)]


def _normalize_key(value: FieldKey) -> FieldKey:
    if isinstance(value, int):
        if value < 0:
            raise ValueError("field row id must be a non-negative integer")
        return value
    key = tuple(value)
    if not key:
        raise ValueError("field row key must not be empty")
    if any(not isinstance(item, int) or item < 0 for item in key):
        raise ValueError("field row key entries must be non-negative integers")
    return key
