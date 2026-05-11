"""Repository for beam profile definitions."""

from __future__ import annotations

from typing import Iterator

from femaster_api.utils import normalize_name

from .profile import Profile


class ProfileRepository:
    """Repository for named beam profile definitions."""

    def __init__(self) -> None:
        self._items: dict[str, Profile] = {}

    def add(self, profile: Profile) -> Profile:
        if not isinstance(profile, Profile):
            raise TypeError("profile must be a Profile")
        self._items[normalize_name(profile.name)] = profile
        return profile

    def get(self, name: str) -> Profile:
        key = normalize_name(name)
        try:
            return self._items[key]
        except KeyError as exc:
            raise KeyError(f"unknown profile: {key}") from exc

    def __getitem__(self, name: str) -> Profile:
        return self.get(name)

    def has(self, value: str | Profile) -> bool:
        if isinstance(value, Profile):
            return any(item is value for item in self._items.values())
        return normalize_name(value) in self._items

    def all(self) -> tuple[Profile, ...]:
        return tuple(self._items[key] for key in sorted(self._items))

    def __iter__(self) -> Iterator[Profile]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._items)
