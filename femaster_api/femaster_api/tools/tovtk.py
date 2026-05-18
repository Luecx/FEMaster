"""VTK export helper for models and result fields."""

from __future__ import annotations

from math import sqrt
from pathlib import Path
from typing import Iterable

from femaster_api.backend import Result
from femaster_api.model import Element, ElementTopology, Field, FieldDomain, Model


def tovtk(
    model: Model,
    result_or_path: Result | str | Path,
    path: str | Path | None = None,
    *,
    loadcase: int | None = None,
    frame: int | None = None,
) -> Path:
    """Write a model and optional result fields to a VTK unstructured grid file."""

    if path is None:
        if isinstance(result_or_path, Result):
            raise TypeError("tovtk(model, result, path) requires an output path")
        result = None
        output_path = Path(result_or_path)
    else:
        if not isinstance(result_or_path, Result):
            raise TypeError("tovtk(model, result, path) expects a Result as second argument")
        result = result_or_path
        output_path = Path(path)
    vtk = _import_vtk()
    grid = vtk.vtkUnstructuredGrid()
    nodes = tuple(model.nodes)
    elements = tuple(model.elements)
    _add_points(vtk, grid, nodes)
    _add_cells(vtk, grid, elements, _node_index_by_id(nodes))

    used = _UsedNames()
    _add_fields(vtk, grid, model.fields, nodes, elements, used, prefix="")
    if result is not None:
        frames = _selected_frames(result, loadcase=loadcase, frame=frame)
        multi_frame = len(frames) > 1
        for loadcase_id, frame_id, fields in frames:
            prefix = f"LC{loadcase_id}_F{frame_id}_" if multi_frame else ""
            _add_fields(vtk, grid, fields.values(), nodes, elements, used, prefix=prefix)

    _write_grid(vtk, grid, output_path)
    return output_path


def _add_points(vtk, grid, nodes) -> None:
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(nodes))
    for index, node in enumerate(nodes):
        points.SetPoint(index, float(node.x), float(node.y), float(node.z))
    grid.SetPoints(points)


def _add_cells(vtk, grid, elements: tuple[Element, ...], node_index: dict[int, int]) -> None:
    for element in elements:
        ids = vtk.vtkIdList()
        for node in element.nodes:
            ids.InsertNextId(node_index[_id(node.id, "node")])
        grid.InsertNextCell(_vtk_cell_type(vtk, element), ids)


def _add_fields(vtk, grid, fields: Iterable[Field], nodes, elements: tuple[Element, ...], used: "_UsedNames", *, prefix: str) -> None:
    node_index = _node_index_by_id(nodes)
    element_index = _element_index_by_id(elements)
    for field in fields:
        if field.domain is FieldDomain.NODE:
            tuples = _indexed_tuples(field, len(nodes), node_index)
            _add_data_array(vtk, grid.GetPointData(), prefix + field.name, tuples, field.cols, used.point)
            _add_vector_derived(vtk, grid.GetPointData(), prefix, field, tuples, used.point)
        elif field.domain is FieldDomain.ELEMENT:
            tuples = _indexed_tuples(field, len(elements), element_index)
            _add_data_array(vtk, grid.GetCellData(), prefix + field.name, tuples, field.cols, used.cell)
            _add_mises_derived(vtk, grid.GetCellData(), prefix, field, tuples, used.cell)
        elif field.domain is FieldDomain.ELEMENT_NODAL:
            point_tuples = _element_nodal_point_average(field, nodes)
            cell_tuples = _element_nodal_cell_tuples(field, elements)
            _add_data_array(vtk, grid.GetCellData(), prefix + field.name + "_ELEMENT_NODAL", cell_tuples, len(cell_tuples[0]), used.cell)
            _add_data_array(vtk, grid.GetPointData(), prefix + field.name, point_tuples, field.cols, used.point)
            _add_vector_derived(vtk, grid.GetPointData(), prefix, field, point_tuples, used.point)
            _add_mises_derived(vtk, grid.GetPointData(), prefix, field, point_tuples, used.point)
        elif field.domain is FieldDomain.ELEMENT_IP:
            cell_tuples = _integration_point_cell_tuples(field, elements)
            _add_data_array(vtk, grid.GetCellData(), prefix + field.name + "_IP", cell_tuples, len(cell_tuples[0]), used.cell)
            averaged = _integration_point_cell_average(field, elements)
            _add_data_array(vtk, grid.GetCellData(), prefix + field.name, averaged, field.cols, used.cell)
            _add_mises_derived(vtk, grid.GetCellData(), prefix, field, averaged, used.cell)
        else:
            tuples = tuple(_clean_values(values, field.cols) for _, values in sorted(field.values.items(), key=lambda item: str(item[0])))
            _add_data_array(vtk, grid.GetFieldData(), prefix + field.name, tuples, field.cols, used.field)


def _indexed_tuples(field: Field, count: int, index_by_id: dict[int, int]) -> tuple[tuple[float, ...], ...]:
    rows = [[0.0] * field.cols for _ in range(count)]
    for key, values in field.values.items():
        if not isinstance(key, int):
            continue
        index = index_by_id.get(key)
        if index is not None:
            rows[index] = list(_clean_values(values, field.cols))
    return tuple(tuple(row) for row in rows)


def _element_nodal_point_average(field: Field, nodes) -> tuple[tuple[float, ...], ...]:
    node_index = _node_index_by_id(nodes)
    sums = [[0.0] * field.cols for _ in nodes]
    counts = [0] * len(nodes)
    for key, values in field.values.items():
        if not isinstance(key, tuple) or len(key) < 2:
            continue
        node_id = key[1]
        index = node_index.get(node_id)
        if index is None:
            continue
        row = _clean_values(values, field.cols)
        for component, value in enumerate(row):
            sums[index][component] += value
        counts[index] += 1
    return tuple(tuple(value / counts[index] if counts[index] else 0.0 for value in row) for index, row in enumerate(sums))


def _element_nodal_cell_tuples(field: Field, elements: tuple[Element, ...]) -> tuple[tuple[float, ...], ...]:
    max_nodes = max((len(element.nodes) for element in elements), default=1)
    width = max_nodes * field.cols
    rows = []
    for element in elements:
        row: list[float] = []
        element_id = _id(element.id, "element")
        for node in element.nodes:
            row.extend(_clean_values(field.values.get((element_id, _id(node.id, "node")), ()), field.cols))
        row.extend([0.0] * (width - len(row)))
        rows.append(tuple(row))
    return tuple(rows) or ((0.0,) * width,)


def _integration_point_cell_tuples(field: Field, elements: tuple[Element, ...]) -> tuple[tuple[float, ...], ...]:
    values_by_element: dict[int, list[tuple[float, ...]]] = {}
    for key, values in field.values.items():
        if isinstance(key, tuple) and key:
            values_by_element.setdefault(key[0], []).append(_clean_values(values, field.cols))
    max_points = max((len(values) for values in values_by_element.values()), default=1)
    width = max_points * field.cols
    rows = []
    for element in elements:
        row: list[float] = []
        for values in values_by_element.get(_id(element.id, "element"), ()):
            row.extend(values)
        row.extend([0.0] * (width - len(row)))
        rows.append(tuple(row))
    return tuple(rows) or ((0.0,) * width,)


def _integration_point_cell_average(field: Field, elements: tuple[Element, ...]) -> tuple[tuple[float, ...], ...]:
    values_by_element: dict[int, list[tuple[float, ...]]] = {}
    for key, values in field.values.items():
        if isinstance(key, tuple) and key:
            values_by_element.setdefault(key[0], []).append(_clean_values(values, field.cols))
    rows = []
    for element in elements:
        values = values_by_element.get(_id(element.id, "element"), ())
        if not values:
            rows.append((0.0,) * field.cols)
            continue
        rows.append(tuple(sum(row[i] for row in values) / len(values) for i in range(field.cols)))
    return tuple(rows)


def _add_vector_derived(vtk, data, prefix: str, field: Field, tuples: tuple[tuple[float, ...], ...], used: set[str]) -> None:
    target = _vector_derived_name(field.name)
    if target is None or field.cols < 3:
        return
    _add_data_array(vtk, data, prefix + target, tuple(row[:3] for row in tuples), 3, used)


def _add_mises_derived(vtk, data, prefix: str, field: Field, tuples: tuple[tuple[float, ...], ...], used: set[str]) -> None:
    if "STRESS" not in field.name.upper() or field.cols < 6:
        return
    rows = tuple((_mises(row),) for row in tuples)
    _add_data_array(vtk, data, prefix + "MISES", rows, 1, used)


def _add_data_array(vtk, data, name: str, tuples: tuple[tuple[float, ...], ...], components: int, used: set[str]) -> str:
    array_name = _unique_name(name, used)
    array = vtk.vtkDoubleArray()
    array.SetName(array_name)
    array.SetNumberOfComponents(components)
    array.SetNumberOfTuples(len(tuples))
    for index, row in enumerate(tuples):
        array.SetTuple(index, _clean_values(row, components))
    data.AddArray(array)
    return array_name


def _selected_frames(result: Result, *, loadcase: int | None, frame: int | None):
    selected = []
    for loadcase_item in result.loadcases:
        if loadcase is not None and loadcase_item.id != loadcase:
            continue
        for frame_item in loadcase_item.frames:
            if frame is not None and frame_item.id != frame:
                continue
            selected.append((loadcase_item.id, frame_item.id, frame_item.fields))
    if not selected:
        raise KeyError("result selection contains no frames")
    return tuple(selected)


def _write_grid(vtk, grid, path: Path) -> None:
    writer = vtk.vtkUnstructuredGridWriter() if path.suffix.lower() == ".vtk" else vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(grid)
    if writer.Write() != 1:
        raise OSError(f"failed to write VTK file: {path}")


def _vtk_cell_type(vtk, element: Element) -> int:
    if element.topology in (ElementTopology.TRUSS2, ElementTopology.BEAM2):
        return vtk.VTK_QUADRATIC_EDGE if len(element.nodes) == 3 else vtk.VTK_LINE
    mapping = {
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
        return mapping[element.topology]
    except KeyError as exc:
        raise ValueError(f"unsupported VTK element topology: {element.topology.value}") from exc


def _vector_derived_name(name: str) -> str | None:
    mapping = {
        "DISPLACEMENT": "DISPLACEMENT_XYZ",
        "MODE_SHAPE": "MODE_SHAPE_XYZ",
        "BUCKLING_MODE": "BUCKLING_MODE_XYZ",
    }
    return mapping.get(name.upper())


def _mises(row: tuple[float, ...]) -> float:
    sxx, syy, szz, sxy, syz, sxz = row[:6]
    return sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2) + 3.0 * (sxy**2 + syz**2 + sxz**2))


def _clean_values(values: Iterable[float | None], cols: int) -> tuple[float, ...]:
    row = tuple(0.0 if value is None else float(value) for value in values)
    return row[:cols] + (0.0,) * max(0, cols - len(row))


def _node_index_by_id(nodes) -> dict[int, int]:
    return {_id(node.id, "node"): index for index, node in enumerate(nodes)}


def _element_index_by_id(elements: tuple[Element, ...]) -> dict[int, int]:
    return {_id(element.id, "element"): index for index, element in enumerate(elements)}


def _id(value: int | None, label: str) -> int:
    if value is None:
        raise ValueError(f"{label} has no repository id")
    return value


def _unique_name(name: str, used: set[str]) -> str:
    if name not in used:
        used.add(name)
        return name
    index = 1
    while f"{name}_{index}" in used:
        index += 1
    unique = f"{name}_{index}"
    used.add(unique)
    return unique


def _import_vtk():
    try:
        import vtk
    except ImportError as exc:
        raise RuntimeError("vtk is required for VTK export") from exc
    return vtk


class _UsedNames:
    def __init__(self) -> None:
        self.point: set[str] = set()
        self.cell: set[str] = set()
        self.field: set[str] = set()
