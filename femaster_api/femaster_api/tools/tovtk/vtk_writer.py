"""Legacy `.vtk` writer based on the Python VTK package."""

from __future__ import annotations

from pathlib import Path

from femaster_api.backend.result_reader import Frame, Result
from femaster_api.model import Model

from .mesh_export import VTK_CELL_TYPES, build_vtk_mesh


class VtkWriter:
    """Write FEMaster model geometry and optional results to legacy VTK."""

    def __init__(
        self,
        model: Model,
        *,
        result: Result | None = None,
        loadcase: str | None = None,
        frame: int = 0,
        binary: bool = False,
    ) -> None:
        self.model    = model
        self.result   = result
        self.loadcase = loadcase
        self.frame    = frame
        self.binary   = binary

    def write(self, path: str | Path) -> None:
        """Write a legacy VTK unstructured-grid file."""

        vtk = _import_vtk()
        grid = _build_grid(vtk, self.model)
        frame = _select_frame(self.result, self.loadcase, self.frame)
        if frame is not None:
            _add_result_fields(vtk, grid, frame)

        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileName(str(path))
        writer.SetInputData(grid)
        if self.binary:
            writer.SetFileTypeToBinary()
        else:
            writer.SetFileTypeToASCII()
        if writer.Write() != 1:
            raise RuntimeError(f"failed to write VTK file: {path}")


def write_vtk(
    model: Model,
    path: str | Path,
    *,
    result: Result | None = None,
    loadcase: str | None = None,
    frame: int = 0,
    binary: bool = False,
) -> None:
    """Write model geometry and optional result fields to a `.vtk` file."""

    VtkWriter(model, result=result, loadcase=loadcase, frame=frame, binary=binary).write(path)


def _build_grid(vtk, model: Model):
    mesh = build_vtk_mesh(model)

    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(mesh.points))
    for index, point in enumerate(mesh.points):
        points.SetPoint(index, *point)

    grid = vtk.vtkUnstructuredGrid()
    grid.SetPoints(points)

    for element in model.elements:
        cell_type = VTK_CELL_TYPES[element.topology]
        ids = vtk.vtkIdList()
        for node in element.nodes:
            ids.InsertNextId(mesh.node_ids[node])
        grid.InsertNextCell(cell_type, ids)

    return grid


def _add_result_fields(vtk, grid, frame: Frame) -> None:
    point_count = grid.GetNumberOfPoints()
    cell_count = grid.GetNumberOfCells()

    for field in frame.fields.values():
        point_rows = _rows_for_count(field.data, point_count)
        if point_rows is not None:
            grid.GetPointData().AddArray(_array_from_rows(vtk, field.name, point_rows))
            continue

        cell_rows = _rows_for_count(field.data, cell_count)
        if cell_rows is not None:
            grid.GetCellData().AddArray(_array_from_rows(vtk, field.name, cell_rows))


def _rows_for_count(
    rows: tuple[tuple[float, ...], ...],
    count: int,
) -> tuple[tuple[float, ...], ...] | None:
    if count == 0 or not rows:
        return None
    if len(rows) == count:
        return rows
    if len(rows) == count + 1:
        return rows[1:]
    return None


def _array_from_rows(vtk, name: str, rows: tuple[tuple[float, ...], ...]):
    component_count = max((len(row) for row in rows), default=1)
    array = vtk.vtkDoubleArray()
    array.SetName(name)
    array.SetNumberOfComponents(component_count)
    array.SetNumberOfTuples(len(rows))

    for row_index, row in enumerate(rows):
        values = tuple(float(value) for value in row)
        padded = values + (0.0,) * (component_count - len(values))
        for component_index, value in enumerate(padded):
            array.SetComponent(row_index, component_index, value)

    return array


def _select_frame(result: Result | None, loadcase: str | None, frame: int) -> Frame | None:
    if result is None or not result.loadcases:
        return None
    if loadcase is None:
        selected_loadcase = next(iter(result.loadcases.values()))
    else:
        selected_loadcase = result.loadcases[loadcase]
    if not selected_loadcase.frames:
        return None
    return selected_loadcase.frames[frame]


def _import_vtk():
    try:
        import vtk
    except ImportError as exc:
        raise RuntimeError("VTK export requires the optional 'vtk' Python package") from exc
    return vtk
