"""Public VTK export API."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

from femaster_api.backend import Result
from femaster_api.model import Model

from .writer_vtk import write_vtk
from .writer_vtkhdf import write_vtkhdf


def tovtk(
    model: Model,
    result_or_path: Result | str | Path,
    path: str | Path | None = None,
    *,
    loadcase: int | None = None,
    frame: int | None = None,
    export_sets: bool = False,
    pack_sets: bool = False,
    export_shell_thickness: bool = False,
    export_couplings: bool = False,
    log_mode: Literal["quiet", "summary", "verbose"] = "quiet",
    log_format: Literal["table", "json"] = "table",
) -> Path:
    result, output_path = _resolve_call(result_or_path, path)
    writer_kwargs = dict(
        export_sets=export_sets,
        pack_sets=pack_sets,
        export_shell_thickness=export_shell_thickness,
        export_couplings=export_couplings,
        log_mode=log_mode,
        log_format=log_format,
    )
    if output_path.suffix.lower() == ".vtkhdf":
        return write_vtkhdf(model, result, output_path, **writer_kwargs)
    return write_vtk(model, result, output_path, loadcase=loadcase, frame=frame, **writer_kwargs)


def _resolve_call(result_or_path: Result | str | Path, path: str | Path | None) -> tuple[Result | None, Path]:
    if path is None:
        if isinstance(result_or_path, Result):
            raise TypeError("tovtk(model, result, path) requires an output path")
        return None, Path(result_or_path)
    if not isinstance(result_or_path, Result):
        raise TypeError("tovtk(model, result, path) expects a Result as second argument")
    return result_or_path, Path(path)
