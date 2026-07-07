"""Small VTK export report renderer."""

from __future__ import annotations

import json


class ExportReport:
    def __init__(self, label: str, *, mode: str, fmt: str) -> None:
        self.label = label
        self.mode = mode
        self.fmt = fmt
        self.rows: list[dict[str, object]] = []

    def add(self, name: str, action: str, target: str, shape: tuple[int, ...] | None, note: str = "") -> None:
        self.rows.append({"field": name, "action": action, "target": target, "shape": shape, "note": note})

    def render(self) -> None:
        if self.mode == "quiet" or not self.rows:
            return
        if self.fmt == "json":
            print(json.dumps({"scope": self.label, "rows": self.rows}, separators=(",", ":")))
            return

        headers = ("Field", "Action", "Target", "Shape", "Note")
        widths = [len(value) for value in headers]
        for row in self.rows:
            values = self._values(row)
            widths = [max(width, len(value)) for width, value in zip(widths, values)]

        print(f"[VTK] {self.label}")
        print("  " + "  ".join(value.ljust(widths[index]) for index, value in enumerate(headers)))
        print("  " + "  ".join("-" * width for width in widths))
        for row in self.rows:
            print("  " + "  ".join(value.ljust(widths[index]) for index, value in enumerate(self._values(row))))

    @staticmethod
    def _values(row: dict[str, object]) -> tuple[str, str, str, str, str]:
        return (
            str(row["field"]),
            str(row["action"]),
            str(row["target"]),
            "-" if row["shape"] is None else str(row["shape"]),
            str(row["note"]),
        )
