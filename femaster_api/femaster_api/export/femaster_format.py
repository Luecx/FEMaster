"""Small formatting helpers for FEMaster deck generation."""

from __future__ import annotations

from collections.abc import Iterable


def number(value: float | int) -> str:
    if isinstance(value, int):
        return str(value)
    return f"{float(value):.12g}"


def csv(values: Iterable[object]) -> str:
    return ", ".join(_token(value) for value in values)


def keyword(command: str, **keys: object) -> str:
    parts = [f"*{command.upper()}"]
    for key, value in keys.items():
        if value is None or value is False:
            continue
        if value is True:
            parts.append(key.upper())
        else:
            parts.append(f"{key.upper()}={value}")
    return ", ".join(parts)


def target_token(value: object) -> str:
    name = getattr(value, "name", None)
    if isinstance(name, str) and name:
        return name
    raw_id = getattr(value, "id", None)
    if raw_id is not None:
        return str(raw_id)
    return str(value)


def id_token(value: object) -> int:
    raw_id = getattr(value, "id", None)
    return int(raw_id if raw_id is not None else value)


def trim_missing(values: Iterable[float | int | None]) -> tuple[float | int | None, ...]:
    result = list(values)
    while result and result[-1] is None:
        result.pop()
    return tuple(result)


def block(lines: Iterable[str]) -> str:
    clean = [line for line in lines if line != ""]
    return "\n".join(clean)


def join_blocks(blocks: Iterable[str]) -> str:
    return "\n\n".join(block for block in blocks if block.strip()) + "\n"


def _token(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float | int):
        return number(value)
    return str(value)
