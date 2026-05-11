"""Repository for material definitions."""

from __future__ import annotations

from typing import Iterator

from femaster_api.utils import normalize_name

from .material import Material


class MaterialRepository:
    """Repository for named material definitions."""

    def __init__(self) -> None:
        self._items: dict[str, Material] = {}

    def add(self, material: Material) -> Material:
        if not isinstance(material, Material):
            raise TypeError("material must be a Material")
        self._items[normalize_name(material.name)] = material
        return material

    def get(self, name: str) -> Material:
        key = normalize_name(name)
        try:
            return self._items[key]
        except KeyError as exc:
            raise KeyError(f"unknown material: {key}") from exc

    def __getitem__(self, name: str) -> Material:
        return self.get(name)

    def has(self, value: str | Material) -> bool:
        if isinstance(value, Material):
            return any(item is value for item in self._items.values())
        return normalize_name(value) in self._items

    def all(self) -> tuple[Material, ...]:
        return tuple(self._items[key] for key in sorted(self._items))

    def __iter__(self) -> Iterator[Material]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._items)
