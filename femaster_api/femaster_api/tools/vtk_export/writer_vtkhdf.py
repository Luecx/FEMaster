"""VTK-HDF time-series writer."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from femaster_api.backend import Frame, LoadCase, Result
from femaster_api.model import Model

from .common import NameRegistry, activate_preferred_arrays, load_vtk
from .fields import write_fields
from .mesh import build_grid
from .report import ExportReport


@dataclass(frozen=True)
class FrameData:
    loadcase: LoadCase | None
    frame: Frame | None
    index: int
    time: float


def write_vtkhdf(
    model: Model,
    result: Result | None,
    out_path: str | Path,
    *,
    export_sets: bool = False,
    pack_sets: bool = False,
    export_shell_thickness: bool = False,
    export_couplings: bool = False,
    log_mode: str = "quiet",
    log_format: str = "table",
) -> Path:
    vtk = load_vtk()
    if not hasattr(vtk, "vtkHDFWriter"):
        raise RuntimeError("this VTK build does not provide vtkHDFWriter")

    base_grid, ctx, names = build_grid(
        vtk,
        model,
        export_sets=export_sets,
        pack_sets=pack_sets,
        export_shell_thickness=export_shell_thickness,
        export_couplings=export_couplings,
    )
    model_report = ExportReport("model fields", mode=log_mode, fmt=log_format)
    write_fields(vtk, base_grid, model.fields, ctx, names, report=model_report)
    model_report.render()

    frames = _frame_data(result)
    algorithm = _temporal_algorithm(vtk, base_grid, ctx, frames, log_mode=log_mode, log_format=log_format)

    out = Path(out_path)
    writer = vtk.vtkHDFWriter()
    writer.SetFileName(str(out))
    writer.SetInputConnection(algorithm.GetOutputPort())
    writer.SetWriteAllTimeSteps(True)
    if writer.Write() != 1:
        raise OSError(f"failed to write VTK-HDF file: {out}")
    return out


def _frame_data(result: Result | None) -> tuple[FrameData, ...]:
    if result is None or not result.loadcases:
        return (FrameData(None, None, 0, 0.0),)
    data: list[FrameData] = []
    counter = 0
    for loadcase in result.loadcases:
        if not loadcase.frames:
            data.append(FrameData(loadcase, None, counter, float(counter)))
            counter += 1
            continue
        for frame in loadcase.frames:
            data.append(FrameData(loadcase, frame, counter, _time_value(frame, counter)))
            counter += 1
    data.sort(key=lambda item: (item.time, item.index))
    for index in range(1, len(data)):
        if data[index].time <= data[index - 1].time:
            data[index] = FrameData(data[index].loadcase, data[index].frame, data[index].index, data[index - 1].time + 1e-12)
    return tuple(data)


def _time_value(frame: Frame, fallback: int) -> float:
    for name in ("TIME", "LAMBDA"):
        field = frame.fields.get(name)
        if field is None:
            continue
        for values in field.values.values():
            if values:
                return float(values[0] or 0.0)
    return float(frame.id if frame.id is not None else fallback)


def _temporal_algorithm(vtk, base_grid, ctx, frames: tuple[FrameData, ...], *, log_mode: str, log_format: str):
    from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase

    class Algorithm(VTKPythonAlgorithmBase):
        def __init__(self):
            super().__init__(nInputPorts=0, nOutputPorts=1, outputType="vtkUnstructuredGrid")
            self.times = [item.time for item in frames]

        def RequestInformation(self, request: Any, in_info: Any, out_info: Any) -> int:
            info = out_info.GetInformationObject(0)
            pipeline = vtk.vtkStreamingDemandDrivenPipeline
            info.Set(pipeline.TIME_STEPS(), self.times, len(self.times))
            info.Set(pipeline.TIME_RANGE(), (self.times[0], self.times[-1]), 2)
            return 1

        def RequestData(self, request: Any, in_info: Any, out_info: Any) -> int:
            info = out_info.GetInformationObject(0)
            pipeline = vtk.vtkStreamingDemandDrivenPipeline
            requested = self.times[0]
            if info.Has(pipeline.UPDATE_TIME_STEP()):
                requested = float(info.Get(pipeline.UPDATE_TIME_STEP()))
            index = min(range(len(self.times)), key=lambda item: abs(self.times[item] - requested))
            frame_data = frames[index]

            grid = vtk.vtkUnstructuredGrid()
            grid.ShallowCopy(base_grid)
            if frame_data.frame is not None:
                report = ExportReport(f"frame={frame_data.frame.id}", mode=log_mode if log_mode == "verbose" else "quiet", fmt=log_format)
                write_fields(vtk, grid, frame_data.frame.fields.values(), ctx, NameRegistry(), report=report)
                report.render()
            activate_preferred_arrays(grid)
            output = vtk.vtkUnstructuredGrid.GetData(out_info, 0)
            output.ShallowCopy(grid)
            return 1

    return Algorithm()
