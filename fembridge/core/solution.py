from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

_NUM_SUFFIX = re.compile(r"^(?P<base>.*?)(?:[_\-]?)(?P<idx>\d+)$", re.IGNORECASE)


def split_series(name: str) -> Tuple[str, Optional[int]]:
    match = _NUM_SUFFIX.match(name)
    if not match:
        return name, None
    base = match.group("base")
    idx_str = match.group("idx")
    try:
        return base, int(idx_str)
    except ValueError:
        return name, None


@dataclass
class Frame:
    index: int
    fields: Dict[str, np.ndarray] = field(default_factory=dict)

    def set_field(self, name: str, data: np.ndarray) -> None:
        self.fields[name] = np.asarray(data)


@dataclass
class StepResult:
    name: str
    loadcase_id: str
    frames: List[Frame] = field(default_factory=list)
    fields: Dict[str, np.ndarray] = field(default_factory=dict)

    def ensure_frame(self, index: int) -> Frame:
        for frame in self.frames:
            if frame.index == index:
                return frame
        frame = Frame(index=index)
        self.frames.append(frame)
        self.frames.sort(key=lambda f: f.index)
        return frame


class Solution:
    """Stores result entries with frames and global fields."""

    def __init__(self) -> None:
        self.steps: List[StepResult] = []
        self.files: Dict[str, Path] = {}

    # ------------------------------------------------------------------
    # Result management helpers (for ASAMI and other engines)
    # ------------------------------------------------------------------
    def clear(self) -> None:
        self.steps.clear()
        self.files.clear()

    def create_result(self, name: str, loadcase_id: str | int) -> StepResult:
        result = StepResult(name=name, loadcase_id=str(loadcase_id))
        self.steps.append(result)
        return result

    def get_result(self, name: str) -> StepResult:
        for result in self.steps:
            if result.name == name:
                return result
        raise KeyError(f"Result '{name}' not found.")

    def ensure_result(self, name: str, loadcase_id: str | int) -> StepResult:
        try:
            return self.get_result(name)
        except KeyError:
            return self.create_result(name, loadcase_id)

    def create_frame(self, result_name: str, frame_index: int) -> Frame:
        result = self.get_result(result_name)
        return result.ensure_frame(frame_index)

    def set_field(
        self,
        result_name: str,
        frame_index: Optional[int],
        field_name: str,
        data: np.ndarray,
        *,
        loadcase_id: Optional[str | int] = None,
    ) -> None:
        if loadcase_id is None:
            result = self.get_result(result_name)
        else:
            result = self.ensure_result(result_name, loadcase_id)

        if frame_index is None:
            result.fields[field_name] = np.asarray(data)
        else:
            frame = result.ensure_frame(frame_index)
            frame.set_field(field_name, data)

    # ------------------------------------------------------------------
    # FEMASTER reader
    # ------------------------------------------------------------------
    def read_res(self, result_path: Path, steps: Optional[Iterable[Any]] = None) -> None:
        self.clear()
        loadcases = self._parse_res_file(result_path)
        step_list = list(steps) if steps is not None else []

        if step_list and len(loadcases) < len(step_list):
            raise ValueError(
                "Result file contains fewer loadcases than model steps. "
                f"Loadcases: {len(loadcases)}, steps: {len(step_list)}"
            )

        for idx, (loadcase_id, fields) in enumerate(loadcases):
            step = step_list[idx] if idx < len(step_list) else None
            name = self._resolve_step_name(step, idx)
            result = self.create_result(name=name, loadcase_id=loadcase_id)
            frame_map: Dict[int, Frame] = {}

            for field_name, data in fields.items():
                base, series_idx = split_series(field_name)
                if series_idx is None:
                    result.fields[field_name] = np.asarray(data)
                else:
                    frame = frame_map.get(series_idx)
                    if frame is None:
                        frame = result.ensure_frame(series_idx)
                        frame_map[series_idx] = frame
                    frame.set_field(base, data)

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------
    def _parse_res_file(self, path: Path) -> List[Tuple[str, Dict[str, np.ndarray]]]:
        loadcases: List[Tuple[str, Dict[str, np.ndarray]]] = []
        current_lc: Optional[str] = None
        current_fields: Dict[str, np.ndarray] = {}
        current_field_name: Optional[str] = None
        matrix_rows: List[List[float]] = []

        with Path(path).open("r", encoding="utf-8") as fh:
            for raw_line in fh:
                line = raw_line.strip()
                if not line:
                    continue

                if line.startswith("LC"):
                    if current_field_name is not None:
                        current_fields[current_field_name] = self._rows_to_array(matrix_rows)
                        current_field_name = None
                        matrix_rows = []
                    if current_lc is not None:
                        loadcases.append((current_lc, current_fields))
                        current_fields = {}
                    current_lc = line.split()[1].strip()

                elif line.startswith("FIELD,"):
                    if current_field_name is not None:
                        current_fields[current_field_name] = self._rows_to_array(matrix_rows)
                    current_field_name = self._extract_field_name(line)
                    matrix_rows = []

                elif line == "END FIELD":
                    if current_field_name is not None:
                        current_fields[current_field_name] = self._rows_to_array(matrix_rows)
                        current_field_name = None
                        matrix_rows = []

                elif current_field_name:
                    matrix_rows.append(self._parse_numbers(line))

            if current_field_name is not None:
                current_fields[current_field_name] = self._rows_to_array(matrix_rows)
            if current_lc is not None:
                loadcases.append((current_lc, current_fields))

        return loadcases

    def _resolve_step_name(self, step: Any, index: int) -> str:
        if step is None:
            return f"STEP_{index + 1}"
        name = getattr(step, "name", None)
        if name:
            return str(name)
        loadcase_type = getattr(step, "loadcase_type", None)
        if loadcase_type:
            return f"{loadcase_type}_{index + 1}"
        return f"{step.__class__.__name__}_{index + 1}"

    @staticmethod
    def _extract_field_name(header_line: str) -> str:
        _, _, name_part = header_line.partition("NAME=")
        name = name_part.split(",", 1)[0].strip()
        if not name:
            raise ValueError(f"Invalid FIELD header: {header_line}")
        return name

    @staticmethod
    def _rows_to_array(rows: List[List[float]]) -> np.ndarray:
        if not rows:
            return np.empty((0, 0), dtype=float)
        return np.array(rows, dtype=float)

    @staticmethod
    def _parse_numbers(line: str) -> List[float]:
        try:
            return [float(chunk) for chunk in line.split()]
        except ValueError as exc:
            raise ValueError(f"Could not parse numeric values from line: '{line}'") from exc

    def __str__(self) -> str:
        if not self.steps:
            return "Solution(results=0)"

        lines: List[str] = ["Solution:"]
        for result in self.steps:
            lines.append(f"- {result.name} (LC {result.loadcase_id})")

            if result.fields:
                lines.append("  Fields:")
                for field_name in sorted(result.fields):
                    data = result.fields[field_name]
                    shape = getattr(data, "shape", None)
                    shape_repr = "(unknown)" if shape is None else str(tuple(shape))
                    lines.append(f"    {field_name}: shape={shape_repr}")

            frames = sorted(result.frames, key=lambda f: f.index)
            if not frames:
                continue

            lines.append("  Frames:")
            total = len(frames)
            if total <= 10:
                leading = frames
                trailing: List[Frame] = []
            else:
                leading = frames[:5]
                trailing = frames[-2:]

            def append_frame(frame: Frame) -> None:
                lines.append(f"    Frame {frame.index}:")
                if not frame.fields:
                    lines.append("      (no fields)")
                    return
                for field_name in sorted(frame.fields):
                    data = frame.fields[field_name]
                    shape = getattr(data, "shape", None)
                    shape_repr = "(unknown)" if shape is None else str(tuple(shape))
                    lines.append(f"      {field_name}: shape={shape_repr}")

            for frame in leading:
                append_frame(frame)

            if trailing:
                omitted = total - (len(leading) + len(trailing))
                plural = "frame" if omitted == 1 else "frames"
                lines.append(f"    ... ({omitted} {plural} omitted)")
                for frame in trailing:
                    append_frame(frame)

        return "\n".join(lines)

    def __repr__(self) -> str:  # pragma: no cover
        return f"Solution(results={len(self.steps)})"
