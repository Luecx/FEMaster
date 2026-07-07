"""Shared VTK export primitives."""

from __future__ import annotations

from dataclasses import dataclass
from math import sqrt
from pathlib import Path
from typing import Iterable, Literal

from femaster_api.model import Element, ElementTopology

VtkTarget = Literal["POINT_DATA", "CELL_DATA", "FIELD_DATA"]


@dataclass(frozen=True)
class MeshContext:
    nodes: tuple[object, ...]
    elements: tuple[Element, ...]
    node_to_point: dict[int, int]
    element_to_cell: dict[int, int]
    cell_count: int

    @property
    def node_count(self) -> int:
        return len(self.nodes)

    @property
    def element_count(self) -> int:
        return len(self.elements)


class NameRegistry:
    def __init__(self) -> None:
        self.point: set[str] = set()
        self.cell: set[str] = set()
        self.field: set[str] = set()

    def for_target(self, target: VtkTarget) -> set[str]:
        if target == "POINT_DATA":
            return self.point
        if target == "CELL_DATA":
            return self.cell
        return self.field


def load_vtk():
    try:
        import vtk
    except ImportError as exc:
        raise RuntimeError("vtk is required for VTK export") from exc
    return vtk


def required_id(value: int | None, label: str) -> int:
    if value is None:
        raise ValueError(f"{label} has no repository id")
    return int(value)


def clean_name(name: str) -> str:
    return str(name).replace("/", "_").replace(" ", "_").replace(".", "_")


def unique_name(name: str, used: set[str]) -> str:
    if name not in used:
        used.add(name)
        return name
    index = 1
    while f"{name}_{index}" in used:
        index += 1
    value = f"{name}_{index}"
    used.add(value)
    return value


def row_values(values: Iterable[float | None], components: int) -> tuple[float, ...]:
    row = tuple(0.0 if value is None else float(value) for value in values)
    return row[:components] + (0.0,) * max(0, components - len(row))


def shape_of(rows: tuple[tuple[float, ...], ...], components: int) -> tuple[int, ...]:
    return (len(rows), components)


def mises(row: tuple[float, ...]) -> float:
    sxx, syy, szz, sxy, syz, sxz = row[:6]
    return sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3.0 * (sxy**2 + syz**2 + sxz**2))


def add_array(vtk, grid, name: str, rows: tuple[tuple[float, ...], ...], target: VtkTarget, names: NameRegistry, *, integer: bool = False) -> str:
    components = len(rows[0]) if rows else 1
    vtk_array = vtk.vtkIntArray() if integer else vtk.vtkDoubleArray()
    array_name = unique_name(clean_name(name), names.for_target(target))
    vtk_array.SetName(array_name)
    vtk_array.SetNumberOfComponents(max(1, components))
    vtk_array.SetNumberOfTuples(len(rows))
    for index, row in enumerate(rows):
        vtk_array.SetTuple(index, row_values(row, components))

    if target == "POINT_DATA":
        grid.GetPointData().AddArray(vtk_array)
    elif target == "CELL_DATA":
        grid.GetCellData().AddArray(vtk_array)
    else:
        grid.GetFieldData().AddArray(vtk_array)
    return array_name


def add_string_array(vtk, grid, name: str, values: list[str], names: NameRegistry) -> None:
    array = vtk.vtkStringArray()
    array.SetName(unique_name(clean_name(name), names.field))
    array.SetNumberOfValues(len(values))
    for index, value in enumerate(values):
        array.SetValue(index, value)
    grid.GetFieldData().AddArray(array)


def activate_preferred_arrays(grid) -> None:
    for data in (grid.GetPointData(), grid.GetCellData()):
        if data is None:
            continue
        for vector_name in ("DISPLACEMENT_XYZ", "MODE_SHAPE_XYZ", "BUCKLING_MODE_XYZ"):
            if data.HasArray(vector_name):
                data.SetActiveVectors(vector_name)
                break
        for scalar_name in ("DISPLACEMENT_MAG", "MODE_SHAPE_MAG", "BUCKLING_MODE_MAG", "MISES", "STRESS_MISES"):
            if data.HasArray(scalar_name):
                data.SetActiveScalars(scalar_name)
                break


def cell_type(vtk, element: Element) -> int:
    if element.topology in (ElementTopology.TRUSS2, ElementTopology.BEAM2):
        return vtk.VTK_QUADRATIC_EDGE if len(element.nodes) == 3 else vtk.VTK_LINE
    types = {
        ElementTopology.TRI3: vtk.VTK_TRIANGLE,
        ElementTopology.QUAD4: vtk.VTK_QUAD,
        ElementTopology.MITC4: vtk.VTK_QUAD,
        ElementTopology.QSPT: vtk.VTK_QUAD,
        ElementTopology.TET4: vtk.VTK_TETRA,
        ElementTopology.PYRAMID5: vtk.VTK_PYRAMID,
        ElementTopology.WEDGE6: vtk.VTK_WEDGE,
        ElementTopology.HEX8: vtk.VTK_HEXAHEDRON,
        ElementTopology.TRI6: vtk.VTK_QUADRATIC_TRIANGLE,
        ElementTopology.QUAD8: vtk.VTK_QUADRATIC_QUAD,
        ElementTopology.TET10: vtk.VTK_QUADRATIC_TETRA,
        ElementTopology.WEDGE15: vtk.VTK_QUADRATIC_WEDGE,
        ElementTopology.HEX20: vtk.VTK_QUADRATIC_HEXAHEDRON,
        ElementTopology.HEX20R: vtk.VTK_QUADRATIC_HEXAHEDRON,
    }
    try:
        return types[element.topology]
    except KeyError as exc:
        raise ValueError(f"unsupported VTK element topology: {element.topology.value}") from exc


def write_unstructured_grid(vtk, grid, path: str | Path) -> None:
    out = Path(path)
    suffix = out.suffix.lower()
    if suffix == ".vtk":
        writer = vtk.vtkUnstructuredGridWriter()
        writer.SetFileTypeToBinary()
    elif suffix == ".vtkhdf":
        if not hasattr(vtk, "vtkHDFWriter"):
            raise RuntimeError("this VTK build does not provide vtkHDFWriter")
        writer = vtk.vtkHDFWriter()
    else:
        writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(out))
    writer.SetInputData(grid)
    if writer.Write() != 1:
        raise OSError(f"failed to write VTK file: {out}")
