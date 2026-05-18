"""Reader for text result files."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

from femaster_api.model.fields import Field, FieldDomain


@dataclass
class Frame:
    id: int = 0
    name: str | None = None
    fields: dict[str, Field] = field(default_factory=dict)

    def field(self, name: str) -> Field:
        return self.fields[name]


@dataclass
class LoadCase:
    id: int
    name: str | None = None
    frames: list[Frame] = field(default_factory=list)

    def frame(self, id: int = 0) -> Frame:
        for frame_item in self.frames:
            if frame_item.id == id:
                return frame_item
        raise KeyError(f"unknown frame id: {id}")

    def field(self, name: str, *, frame: int = 0) -> Field:
        return self.frame(frame).field(name)


@dataclass
class Result:
    loadcases: list[LoadCase] = field(default_factory=list)

    def loadcase(self, id: int = 1) -> LoadCase:
        for loadcase in self.loadcases:
            if loadcase.id == id:
                return loadcase
        raise KeyError(f"unknown loadcase id: {id}")

    def field(self, name: str, *, loadcase: int = 1, frame: int = 0) -> Field:
        return self.loadcase(loadcase).field(name, frame=frame)


class ResultReader:
    """Parse `.res` text into loadcase, frame, and field objects."""

    def read(self, path: str | Path) -> Result:
        return self.parse(Path(path).read_text(encoding="utf-8"))

    def parse(self, text: str) -> Result:
        result = Result()
        current_loadcase: LoadCase | None = None
        current_frame: Frame | None = None
        lines = iter(enumerate(text.splitlines(), start=1))

        for _, raw in lines:
            line = raw.strip()
            if not line or _is_comment(line):
                continue
            upper = line.upper()
            if upper.startswith(("LC ", "LOADCASE ")):
                parts = line.split(maxsplit=2)
                current_loadcase = LoadCase(int(parts[1]), parts[2] if len(parts) > 2 else None)
                current_frame = None
                result.loadcases.append(current_loadcase)
                continue
            if upper.startswith("FRAME"):
                if current_loadcase is None:
                    current_loadcase = LoadCase(1, None)
                    result.loadcases.append(current_loadcase)
                current_frame = _parse_frame(line)
                current_loadcase.frames.append(current_frame)
                continue
            if upper.startswith("FIELD"):
                if current_loadcase is None:
                    current_loadcase = LoadCase(1, None)
                    result.loadcases.append(current_loadcase)
                if current_frame is None:
                    current_frame = Frame()
                    current_loadcase.frames.append(current_frame)
                field_item = _read_field(line, lines)
                current_frame.fields[field_item.name] = field_item

        return result


def _parse_frame(header: str) -> Frame:
    parts = header.split(maxsplit=2)
    if len(parts) == 1:
        return Frame()
    try:
        frame_id = int(parts[1])
        name = parts[2] if len(parts) > 2 else None
    except ValueError:
        frame_id = 0
        name = header.partition(" ")[2] or None
    return Frame(frame_id, name)


def _read_field(header: str, lines: Iterator[tuple[int, str]]) -> Field:
    command, keys = _parse_field_header(header)
    rows = int(keys["ROWS"])
    name = keys.get("NAME", command)
    cols = int(keys["COLS"]) if "COLS" in keys else None
    index_cols = int(keys.get("INDEX_COLS", 0))
    value_cols = int(keys["VALUE_COLS"]) if "VALUE_COLS" in keys else None
    total_cols = cols if cols is not None else index_cols + (value_cols or 0)
    if value_cols is None:
        value_cols = total_cols - index_cols if total_cols is not None else 0

    field_item = Field(
        name=name,
        domain=_field_domain(keys, index_cols),
        cols=value_cols,
        fill="NONE",
    )
    read_rows = 0
    for _, raw in lines:
        line = raw.strip()
        if not line or _is_comment(line):
            continue
        if line.upper().startswith("END FIELD"):
            break
        values = tuple(float(token) for token in line.replace(",", " ").split())
        if total_cols and len(values) != total_cols:
            raise ValueError(f"field {name!r} expected {total_cols} columns, got {len(values)}")
        key = _field_key(values, index_cols, read_rows)
        field_item = field_item.set(key, values[index_cols:])
        read_rows += 1
        if read_rows == rows:
            break

    return field_item


def _field_key(values: tuple[float, ...], index_cols: int, fallback: int) -> int | tuple[int, ...]:
    if index_cols == 0:
        return fallback
    index = tuple(int(value) for value in values[:index_cols])
    return index[0] if len(index) == 1 else index


def _field_domain(keys: dict[str, str], index_cols: int) -> FieldDomain:
    raw = keys.get("TYPE") or keys.get("DOMAIN")
    if raw:
        normalized = raw.replace("_", "").upper()
        if normalized == "IP":
            return FieldDomain.ELEMENT_IP
        for domain in FieldDomain:
            if domain.value.replace("_", "").upper() == normalized or domain.name.replace("_", "").upper() == normalized:
                return domain
    if index_cols > 0:
        return FieldDomain.ELEMENT_NODAL
    return FieldDomain.UNKNOWN


def _parse_field_header(header: str) -> tuple[str, dict[str, str]]:
    normalized = header.replace(",", " ")
    parts = normalized.split()
    if len(parts) >= 2 and "=" not in parts[1]:
        command = parts[1]
        tokens = parts[2:]
    else:
        command = "FIELD"
        tokens = parts[1:]
    keys: dict[str, str] = {}
    for token in tokens:
        if "=" in token:
            key, value = token.split("=", 1)
            keys[key.upper()] = value
    return command, keys


def _is_comment(line: str) -> bool:
    return line.startswith(("#", "!", "//", "**"))
