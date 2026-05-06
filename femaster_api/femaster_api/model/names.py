"""Name validation helpers for model-owned named objects."""

from __future__ import annotations


def require_name(value: str, *, label: str = "name") -> str:
    if not isinstance(value, str):
        raise TypeError(f"{label} must be a string")
    stripped = value.strip()
    if not stripped:
        raise ValueError(f"{label} must not be empty")
    if any(ch.isspace() for ch in stripped):
        raise ValueError(f"{label} must not contain whitespace: {value!r}")
    return stripped
