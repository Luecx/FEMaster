"""Name normalization helpers."""

from __future__ import annotations


def normalize_name(value: str, *, label: str = "name") -> str:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string")
    name = value.strip()
    if not name:
        raise ValueError(f"{label} must not be empty")
    if any(ch.isspace() for ch in name):
        raise ValueError(f"{label} must not contain whitespace: {value!r}")
    return name
