"""Small formatting helpers for FEMaster keyword input.

Formatting stays intentionally explicit. Model objects construct their own
blocks and use these helpers only for consistent keyword, CSV and block syntax.
"""

from __future__ import annotations

from collections.abc import Iterable


def keyword(name: str, **keys: object) -> str:
    """Return one FEMaster keyword line and omit keys whose value is None."""

    parts = [f"*{name.upper()}"]
    for key, value in keys.items():
        if value is None:
            continue
        if isinstance(value, bool):
            value = "YES" if value else "NO"
        parts.append(f"{key.upper()}={value}")
    return ", ".join(parts)


def csv(values: Iterable[object]) -> str:
    """Return one comma-separated FEMaster data line."""

    return ", ".join("" if value is None else str(value) for value in values)


def block(lines: Iterable[str]) -> str:
    """Join non-empty lines into one keyword block."""

    return "\n".join(line for line in lines if line)


def blocks(items: Iterable[str]) -> str:
    """Join non-empty keyword blocks with one separating blank line."""

    return "\n\n".join(item for item in items if item)
