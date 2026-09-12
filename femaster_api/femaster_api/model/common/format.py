"""Small formatting helpers for native FEMaster keyword decks."""

from __future__ import annotations

from collections.abc import Iterable


def format_value(value: object) -> str:
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, float):
        return repr(value)
    return str(value)


def csv(values: Iterable[object]) -> str:
    return ", ".join("" if value is None else format_value(value) for value in values)


def keyword(name: str, **keys: object) -> str:
    result = f"*{name}"
    for key, value in keys.items():
        if value is None:
            continue
        result += f", {key}={format_value(value)}"
    return result


def block(lines: Iterable[str]) -> str:
    return "\n".join(line for line in lines if line != "")


def blocks(items: Iterable[str]) -> str:
    return "\n\n".join(item for item in items if item)
