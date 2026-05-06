"""Repository for generic fields."""

from __future__ import annotations

from typing import Iterator

from femaster_api.model.names import require_name

from .field import Field


class FieldRepository:
    """Repository for named generic fields."""

    def __init__(self) -> None:
        self._items: dict[str, Field] = {}

    def add(self, field: Field) -> Field:
        if not isinstance(field, Field):
            raise TypeError("field must be a Field")
        if field.cols <= 0:
            raise ValueError("field cols must be positive")
        if field.cols > 64:
            raise ValueError("FEMaster fields are limited to 64 columns")
        self._items[require_name(field.name)] = field
        return field

    def get(self, name: str) -> Field:
        key = require_name(name)
        try:
            return self._items[key]
        except KeyError as exc:
            raise KeyError(f"unknown field: {key}") from exc

    def __getitem__(self, name: str) -> Field:
        return self.get(name)

    def has(self, value: str | Field) -> bool:
        if isinstance(value, Field):
            return any(item is value for item in self._items.values())
        return require_name(value) in self._items

    def all(self) -> tuple[Field, ...]:
        return tuple(self._items[key] for key in sorted(self._items))

    def __iter__(self) -> Iterator[Field]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._items)
