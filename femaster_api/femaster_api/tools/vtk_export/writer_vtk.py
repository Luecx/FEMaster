"""Static VTK/VTP-style unstructured-grid writer."""

from __future__ import annotations

from pathlib import Path

from femaster_api.backend import Result
from femaster_api.model import Model

from .common import activate_preferred_arrays, load_vtk, write_unstructured_grid
from .fields import write_fields
from .mesh import build_grid
from .report import ExportReport


def write_vtk(
    model: Model,
    result: Result | None,
    out_path: str | Path,
    *,
    loadcase: int | None = None,
    frame: int | None = None,
    export_sets: bool = False,
    pack_sets: bool = False,
    export_shell_thickness: bool = False,
    export_couplings: bool = False,
    log_mode: str = "quiet",
    log_format: str = "table",
) -> Path:
    vtk = load_vtk()
    grid, ctx, names = build_grid(
        vtk,
        model,
        export_sets=export_sets,
        pack_sets=pack_sets,
        export_shell_thickness=export_shell_thickness,
        export_couplings=export_couplings,
    )

    model_report = ExportReport("model fields", mode=log_mode, fmt=log_format)
    write_fields(vtk, grid, model.fields, ctx, names, report=model_report)
    model_report.render()

    if result is not None:
        selected = _selected_frames(result, loadcase=loadcase, frame=frame)
        multi = len(selected) > 1
        for loadcase_id, frame_item in selected:
            prefix = f"LC{loadcase_id}_F{frame_item.id}_" if multi else ""
            report = ExportReport(f"loadcase={loadcase_id} frame={frame_item.id}", mode=log_mode, fmt=log_format)
            write_fields(vtk, grid, frame_item.fields.values(), ctx, names, report=report, prefix=prefix)
            report.render()

    activate_preferred_arrays(grid)
    write_unstructured_grid(vtk, grid, out_path)
    return Path(out_path)


def _selected_frames(result: Result, *, loadcase: int | None, frame: int | None):
    selected = []
    for loadcase_item in result.loadcases:
        if loadcase is not None and loadcase_item.id != loadcase:
            continue
        for frame_item in loadcase_item.frames:
            if frame is not None and frame_item.id != frame:
                continue
            selected.append((loadcase_item.id, frame_item))
    if not selected:
        raise KeyError("result selection contains no frames")
    return tuple(selected)
