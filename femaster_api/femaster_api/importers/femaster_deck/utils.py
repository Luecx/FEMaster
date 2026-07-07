"""Tokenizer utilities for FEMaster input decks."""

from __future__ import annotations

from dataclasses import dataclass
import csv
from pathlib import Path
from typing import Iterable, Mapping, TextIO


class FEMasterInputError(RuntimeError):
    """Raised for malformed or unsupported FEMaster input constructs."""


@dataclass(frozen=True)
class Header:
    keyword: str
    params: Mapping[str, str]
    raw: str


class LineStream:
    """Cursor over non-empty, non-comment input deck lines."""

    def __init__(self, lines: Iterable[str]) -> None:
        self._lines = iter(lines)
        self._buffer: str | None = None

    def peek(self) -> str | None:
        if self._buffer is None:
            self._buffer = self._next_line()
        return self._buffer

    def pop(self) -> str | None:
        line = self.peek()
        self._buffer = None
        return line

    def _next_line(self) -> str | None:
        for raw in self._lines:
            line = strip_inline_comment(raw).strip()
            if not line or line.startswith("**"):
                continue
            return line
        return None


def read_lines(source: str | Path | TextIO | Iterable[str]) -> list[str]:
    if isinstance(source, (str, Path)):
        return Path(source).read_text(encoding="utf-8").splitlines()
    if hasattr(source, "read"):
        return source.readlines()
    return list(source)


def strip_inline_comment(line: str) -> str:
    marker = line.find("**")
    return line if marker < 0 else line[:marker]


def parse_header(line: str) -> Header:
    if not line.startswith("*"):
        raise FEMasterInputError(f"expected command header, got {line!r}")
    parts = next(csv.reader([line[1:]], skipinitialspace=True))
    if not parts:
        raise FEMasterInputError("empty command header")
    keyword = parts[0].strip().replace(" ", "").upper()
    params: dict[str, str] = {}
    for token in parts[1:]:
        token = token.strip()
        if not token:
            continue
        if "=" in token:
            key, value = token.split("=", 1)
            params[key.strip().replace(" ", "_").upper()] = value.strip()
        else:
            params[token.replace(" ", "_").upper()] = "TRUE"
    return Header(keyword, params, line)


def parse_csv(line: str) -> list[str]:
    return [token.strip() for token in next(csv.reader([line], skipinitialspace=True))]


def truthy(value: str | None) -> bool:
    return value is not None and value.strip().upper() in {"", "1", "TRUE", "YES", "ON"}


def falsey(value: str | None) -> bool:
    return value is not None and value.strip().upper() in {"0", "FALSE", "NO", "OFF"}


def numbers(line: str) -> list[float]:
    return [float(token) for token in parse_csv(line) if token]


def side_token(value: str) -> str:
    token = value.strip().upper()
    if token.startswith("S"):
        token = token[1:]
    try:
        return f"S{int(token)}"
    except ValueError as exc:
        raise FEMasterInputError(f"invalid surface side identifier {value!r}") from exc


def normalize_loadcase_type(value: str | None) -> str:
    if not value:
        raise FEMasterInputError("*LOADCASE requires TYPE")
    return value.replace("_", "").replace(" ", "").upper()
