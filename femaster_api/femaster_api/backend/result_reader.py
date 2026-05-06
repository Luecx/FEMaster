"""Reader for FEMaster text result files."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True, slots=True)
class ResultField:
    name: str
    rows: int
    cols: int | None = None
    index_cols: int = 0
    value_cols: int | None = None
    data: tuple[tuple[float, ...], ...] = ()

    def row(self, index: int) -> tuple[float, ...]:
        return self.data[index]


@dataclass
class LoadCaseResult:
    id: int
    name: str | None = None
    fields: dict[str, ResultField] = field(default_factory=dict)

    def field(self, name: str) -> ResultField:
        return self.fields[name]


@dataclass
class FEMasterResult:
    loadcases: list[LoadCaseResult] = field(default_factory=list)

    def loadcase(self, id: int = 1) -> LoadCaseResult:
        for loadcase in self.loadcases:
            if loadcase.id == id:
                return loadcase
        raise KeyError(f"unknown loadcase id: {id}")

    def field(self, name: str, *, loadcase: int = 1) -> ResultField:
        return self.loadcase(loadcase).field(name)


class ResultReader:
    """Parse FEMaster `.res` files into loadcase and field objects."""

    def read(self, path: str | Path) -> FEMasterResult:
        return self.parse(Path(path).read_text(encoding="utf-8"))

    def parse(self, text: str) -> FEMasterResult:
        result = FEMasterResult()
        current_lc: LoadCaseResult | None = None
        lines = iter(enumerate(text.splitlines(), start=1))

        for _, raw in lines:
            line = raw.strip()
            if not line or _is_comment(line):
                continue
            upper = line.upper()
            if upper.startswith("LC "):
                parts = line.split(maxsplit=2)
                current_lc = LoadCaseResult(int(parts[1]), parts[2] if len(parts) > 2 else None)
                result.loadcases.append(current_lc)
                continue
            if upper.startswith("LOADCASE "):
                parts = line.split(maxsplit=2)
                current_lc = LoadCaseResult(int(parts[1]), parts[2] if len(parts) > 2 else None)
                result.loadcases.append(current_lc)
                continue
            if upper.startswith("FIELD"):
                if current_lc is None:
                    current_lc = LoadCaseResult(1, None)
                    result.loadcases.append(current_lc)
                field = _read_field(line, lines)
                current_lc.fields[field.name] = field

        return result


def _read_field(header: str, lines) -> ResultField:
    command, keys = _parse_field_header(header)
    rows = int(keys["ROWS"])
    name = keys.get("NAME", command)
    cols = int(keys["COLS"]) if "COLS" in keys else None
    index_cols = int(keys.get("INDEX_COLS", 0))
    value_cols = int(keys["VALUE_COLS"]) if "VALUE_COLS" in keys else None
    total_cols = cols if cols is not None else index_cols + (value_cols or 0)

    data: list[tuple[float, ...]] = []
    for _, raw in lines:
        line = raw.strip()
        if not line or _is_comment(line):
            continue
        if line.upper().startswith("END FIELD"):
            break
        values = tuple(float(token) for token in line.replace(",", " ").split())
        if total_cols and len(values) != total_cols:
            raise ValueError(f"field {name!r} expected {total_cols} columns, got {len(values)}")
        data.append(values)
        if len(data) == rows:
            continue

    return ResultField(name, rows, cols, index_cols, value_cols, tuple(data))


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
