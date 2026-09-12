"""Importer for native FEMaster RES text results."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

from ..model import Field, FieldDomain, FieldType, Frame, LoadCase, Result


class ResImporter:
    """Import native FEMaster RES result files."""

    def import_file(self, path: str | Path) -> Result:
        return self.import_text(Path(path).read_text(encoding="utf-8"))

    def import_text(self, text: str) -> Result:
        result = Result()
        current_loadcase: LoadCase | None = None
        current_frame: Frame | None = None
        lines = iter(enumerate(text.splitlines(), start=1))

        for _, raw in lines:
            line = raw.strip()
            if not line or self._is_comment(line):
                continue

            upper = line.upper()

            if upper.startswith(("LC ", "LOADCASE ")):
                parts = line.split(maxsplit=2)
                current_loadcase = LoadCase(
                    int(parts[1]),
                    parts[2] if len(parts) > 2 else None,
                )
                current_frame = None
                result.add(current_loadcase)
                continue

            if upper.startswith("FRAME"):
                if current_loadcase is None:
                    current_loadcase = result.add(LoadCase(1))
                current_frame = self._parse_frame(line)
                current_loadcase.add(current_frame)
                continue

            if upper.startswith("FIELD"):
                if current_loadcase is None:
                    current_loadcase = result.add(LoadCase(1))
                if current_frame is None:
                    current_frame = current_loadcase.add(Frame())
                current_frame.add(self._read_field(line, lines))

        return result

    @staticmethod
    def _parse_frame(header: str) -> Frame:
        parts = header.split(maxsplit=2)
        if len(parts) == 1:
            return Frame()
        try:
            return Frame(int(parts[1]), parts[2] if len(parts) > 2 else None)
        except ValueError:
            return Frame(0, header.partition(" ")[2] or None)

    def _read_field(
        self,
        header: str,
        lines: Iterator[tuple[int, str]],
    ) -> Field:
        command, keys = self._parse_field_header(header)
        rows = int(keys.get("ROWS", "0"))
        name = keys.get("NAME", command)
        index_cols = int(keys.get("INDEX_COLS", "0"))

        if "VALUE_COLS" in keys:
            value_cols = int(keys["VALUE_COLS"])
        elif "COLS" in keys:
            value_cols = int(keys["COLS"])
        else:
            value_cols = 0

        components = tuple(
            item.strip()
            for item in keys.get("COMPONENTS", "").split(";")
            if item.strip()
        )
        if not components:
            components = tuple(
                f"C{index + 1}" for index in range(value_cols)
            )

        field = Field(
            name,
            self._field_domain(keys, index_cols),
            components,
            type=FieldType.from_name(name),
        )

        read_rows = 0
        for _, raw in lines:
            line = raw.strip()
            if not line or self._is_comment(line):
                continue
            if line.upper().startswith("END FIELD"):
                break

            tokens = line.replace(",", " ").split()
            expected = index_cols + value_cols
            if len(tokens) < expected:
                raise ValueError(
                    f"field {name!r} expected at least {expected} columns, "
                    f"got {len(tokens)}"
                )

            key = self._field_key(tokens[:index_cols], read_rows)
            values = tuple(
                float(token)
                for token in tokens[index_cols:index_cols + value_cols]
            )
            field.set(key, values)
            read_rows += 1

            if rows and read_rows >= rows:
                break

        if rows and read_rows != rows:
            raise ValueError(
                f"field {name!r} expected {rows} rows, got {read_rows}"
            )

        return field

    @classmethod
    def _field_key(
        cls,
        tokens: list[str],
        fallback: int,
    ):
        if not tokens:
            return fallback

        values = tuple(cls._semantic_id(token) for token in tokens)
        return values[0] if len(values) == 1 else values

    @staticmethod
    def _semantic_id(token: str):
        try:
            return int(token)
        except ValueError:
            return token

    @staticmethod
    def _field_domain(
        keys: dict[str, str],
        index_cols: int,
    ) -> FieldDomain:
        raw = keys.get("TYPE") or keys.get("DOMAIN")
        if raw:
            normalized = raw.replace("_", "").upper()
            aliases = {
                "IP": FieldDomain.ELEMENT_IP,
                "MP": FieldDomain.ELEMENT_MP,
            }
            if normalized in aliases:
                return aliases[normalized]

            for domain in FieldDomain:
                if normalized in {
                    domain.value.replace("_", "").upper(),
                    domain.name.replace("_", "").upper(),
                }:
                    return domain

        if index_cols > 0:
            return FieldDomain.ELEMENT_NODAL
        return FieldDomain.UNKNOWN

    @staticmethod
    def _parse_field_header(
        header: str,
    ) -> tuple[str, dict[str, str]]:
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

    @staticmethod
    def _is_comment(line: str) -> bool:
        return line.startswith(("#", "!", "//", "**"))
