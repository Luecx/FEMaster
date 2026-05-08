"""Reader for FEMaster text result files."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True, slots=True)
class ResultField:
    """Tabular result field stored in one result frame."""

    name: str
    rows: int
    cols: int | None = None
    index_cols: int = 0
    value_cols: int | None = None
    data: tuple[tuple[float, ...], ...] = ()

    def row(self, index: int) -> tuple[float, ...]:
        return self.data[index]


@dataclass
class Frame:
    """One result frame from a static, transient, modal, or buckling loadcase."""

    fields: dict[str, ResultField] = field(default_factory=dict)
    time: float | None = None
    mode: int | None = None
    eigenvalue: float | None = None
    frequency: float | None = None
    load_factor: float | None = None

    def field(self, name: str) -> ResultField:
        return self.fields[name]


@dataclass
class LoadcaseResult:
    """Result data for one named loadcase."""

    name: str
    id: int | None = None
    frames: list[Frame] = field(default_factory=list)

    @property
    def fields(self) -> dict[str, ResultField]:
        if not self.frames:
            self.frames.append(Frame())
        return self.frames[0].fields

    def field(self, name: str, *, frame: int = 0) -> ResultField:
        return self.frames[frame].fields[name]


LoadCaseResult = LoadcaseResult


@dataclass
class Result:
    """Top-level FEMaster result object."""

    loadcases: dict[str, LoadcaseResult] = field(default_factory=dict)

    def loadcase(self, name: str | int = 1) -> LoadcaseResult:
        if isinstance(name, int):
            for loadcase in self.loadcases.values():
                if loadcase.id == name:
                    return loadcase
            raise KeyError(f"unknown loadcase id: {name}")
        return self.loadcases[name]

    def field(self, name: str, *, loadcase: str | int = 1, frame: int = 0) -> ResultField:
        return self.loadcase(loadcase).field(name, frame=frame)


FEMasterResult = Result


@dataclass
class _ResultCursor:
    result: Result
    current_loadcase: LoadcaseResult | None = None
    current_frame: Frame | None = None

    def use_loadcase(self, id: int | None, name: str | None) -> None:
        key = name or (f"LC{id}" if id is not None else "LC1")
        loadcase = LoadcaseResult(key, id)
        loadcase.frames.append(Frame())
        self.result.loadcases[key] = loadcase
        self.current_loadcase = loadcase
        self.current_frame = loadcase.frames[0]

    def use_frame(self, **metadata: float | int | None) -> None:
        if self.current_loadcase is None:
            self.use_loadcase(1, None)
        frame = Frame(**metadata)
        if len(self.current_loadcase.frames) == 1 and not self.current_loadcase.frames[0].fields:
            self.current_loadcase.frames[0] = frame
            self.current_frame = frame
            return
        self.current_loadcase.frames.append(frame)
        self.current_frame = frame

    def add_field(self, field: ResultField) -> None:
        if self.current_loadcase is None:
            self.use_loadcase(1, None)
        if self.current_frame is None:
            self.current_frame = Frame()
            self.current_loadcase.frames.append(self.current_frame)
        self.current_frame.fields[field.name] = field


class ResultReader:
    """Parse FEMaster `.res` files into loadcase, frame, and field objects."""

    def read(self, path: str | Path) -> Result:
        return self.parse(Path(path).read_text(encoding="utf-8"))

    def parse(self, text: str) -> Result:
        cursor = _ResultCursor(Result())
        lines = iter(enumerate(text.splitlines(), start=1))

        for _, raw in lines:
            line = raw.strip()
            if not line or _is_comment(line):
                continue
            upper = line.upper()
            if upper.startswith("LC "):
                parts = line.split(maxsplit=2)
                cursor.use_loadcase(int(parts[1]), parts[2] if len(parts) > 2 else None)
                continue
            if upper.startswith("LOADCASE "):
                parts = line.split(maxsplit=2)
                cursor.use_loadcase(int(parts[1]), parts[2] if len(parts) > 2 else None)
                continue
            if upper.startswith("FRAME"):
                cursor.use_frame(**_parse_frame_header(line))
                continue
            if upper.startswith("FIELD"):
                cursor.add_field(_read_field(line, lines))

        return cursor.result


def _parse_frame_header(header: str) -> dict[str, float | int | None]:
    keys = _parse_key_values(header)
    metadata: dict[str, float | int | None] = {}
    if "TIME" in keys:
        metadata["time"] = float(keys["TIME"])
    if "MODE" in keys:
        metadata["mode"] = int(keys["MODE"])
    if "EIGENVALUE" in keys:
        metadata["eigenvalue"] = float(keys["EIGENVALUE"])
    if "FREQUENCY" in keys:
        metadata["frequency"] = float(keys["FREQUENCY"])
    if "LOAD_FACTOR" in keys:
        metadata["load_factor"] = float(keys["LOAD_FACTOR"])
    if "LOADFACTOR" in keys:
        metadata["load_factor"] = float(keys["LOADFACTOR"])
    return metadata


def _parse_key_values(line: str) -> dict[str, str]:
    keys: dict[str, str] = {}
    for token in line.replace(",", " ").split()[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            keys[key.upper()] = value
    return keys


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
