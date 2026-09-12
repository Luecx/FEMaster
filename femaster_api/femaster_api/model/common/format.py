"""Low-level text-format helpers for native FEMaster input export.

This module contains only formatting functions.  It deliberately has no model
classes and no knowledge about FEM semantics.  Exportable model objects decide
which keywords and rows they need; these helpers only turn those decisions into
deterministic deck text.

Keeping the formatting functions here avoids a serializer framework while still
ensuring identical comma handling, optional keyword handling and block joining
throughout the public model API.
"""

from __future__ import annotations

from collections.abc import Iterable


def scalar(value: object) -> str:
    """Convert one Python value to the textual form used in deck rows."""

    if value is None:
        return ""
    if isinstance(value, bool):
        return "1" if value else "0"
    return str(value)


def csv(values: Iterable[object]) -> str:
    """Join one FEMaster data row without altering semantic values."""

    return ", ".join(scalar(value) for value in values)


def keyword(name: str, **keys: object) -> str:
    """Build one ``*KEYWORD`` line from non-``None`` keyword parameters."""

    result = f"*{name}"
    rendered = [
        f"{key}={scalar(value)}"
        for key, value in keys.items()
        if value is not None
    ]
    if rendered:
        result += ", " + ", ".join(rendered)
    return result


def block(lines: Iterable[str]) -> str:
    """Join non-empty logical lines into one keyword block."""

    return "\n".join(line for line in lines if line)


def blocks(items: Iterable[str]) -> str:
    """Join non-empty keyword blocks with one visual separator line."""

    return "\n\n".join(item for item in items if item)
