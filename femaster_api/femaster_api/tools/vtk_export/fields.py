"""Field routing for VTK exports."""

from __future__ import annotations

from math import sqrt
from typing import Iterable

from femaster_api.model import Field, FieldDomain

from .common import MeshContext, NameRegistry, VtkTarget, add_array, clean_name, mises, required_id, row_values, shape_of
from .report import ExportReport


def write_fields(vtk, grid, fields: Iterable[Field], ctx: MeshContext, names: NameRegistry, *, report: ExportReport | None = None, prefix: str = "") -> None:
    for field in fields:
        exported = write_field(vtk, grid, field, ctx, names, report=report, prefix=prefix)
        for target, array_name, rows, components in exported:
            add_derived(vtk, grid, field.name, array_name, rows, components, target, names, report)


def write_field(vtk, grid, field: Field, ctx: MeshContext, names: NameRegistry, *, report: ExportReport | None, prefix: str = "") -> list[tuple[VtkTarget, str, tuple[tuple[float, ...], ...], int]]:
    name = clean_name(prefix + field.name)
    if field.domain is FieldDomain.NODE:
        rows = _indexed_rows(field, ctx.node_count, ctx.node_to_point)
        stored = add_array(vtk, grid, name, rows, "POINT_DATA", names)
        _report(report, stored, "WRITE", "POINT_DATA", shape_of(rows, field.cols))
        return [("POINT_DATA", stored, rows, field.cols)]
    if field.domain is FieldDomain.ELEMENT:
        rows = _indexed_rows(field, ctx.cell_count, ctx.element_to_cell)
        stored = add_array(vtk, grid, name, rows, "CELL_DATA", names)
        _report(report, stored, "WRITE", "CELL_DATA", shape_of(rows, field.cols), "element id mapped to cell id")
        return [("CELL_DATA", stored, rows, field.cols)]
    if field.domain is FieldDomain.ELEMENT_NODAL:
        point_rows = _element_nodal_average(field, ctx)
        cell_rows = _element_nodal_packed(field, ctx)
        point_name = add_array(vtk, grid, name, point_rows, "POINT_DATA", names)
        cell_name = add_array(vtk, grid, name + "_ELEMENT_NODAL", cell_rows, "CELL_DATA", names)
        _report(report, point_name, "WRITE", "POINT_DATA", shape_of(point_rows, field.cols), "element-nodal average")
        _report(report, cell_name, "WRITE", "CELL_DATA", shape_of(cell_rows, len(cell_rows[0]) if cell_rows else 0), "packed element-node values")
        return [("POINT_DATA", point_name, point_rows, field.cols)]
    if field.domain is FieldDomain.ELEMENT_IP:
        average = _integration_point_average(field, ctx)
        packed = _integration_point_packed(field, ctx)
        avg_name = add_array(vtk, grid, name, average, "CELL_DATA", names)
        pack_name = add_array(vtk, grid, name + "_IP", packed, "CELL_DATA", names)
        _report(report, avg_name, "WRITE", "CELL_DATA", shape_of(average, field.cols), "integration-point average")
        _report(report, pack_name, "WRITE", "CELL_DATA", shape_of(packed, len(packed[0]) if packed else 0), "packed integration-point values")
        return [("CELL_DATA", avg_name, average, field.cols)]

    routed = _route_unknown(field, ctx)
    if routed is None:
        _report(report, name, "IGNORE", "-", None, "unknown domain could not be matched")
        return []
    target, rows = routed
    stored = add_array(vtk, grid, name, rows, target, names)
    _report(report, stored, "WRITE", target, shape_of(rows, field.cols), "inferred from field keys")
    return [(target, stored, rows, field.cols)]


def add_derived(vtk, grid, field_name: str, array_name: str, rows: tuple[tuple[float, ...], ...], components: int, target: VtkTarget, names: NameRegistry, report: ExportReport | None) -> None:
    if target == "FIELD_DATA":
        return
    upper = field_name.upper()
    if upper in {"DISPLACEMENT", "MODE_SHAPE", "BUCKLING_MODE", "EXTERNAL_FORCES"} and components >= 3:
        vectors = tuple(row[:3] for row in rows)
        vector_name = add_array(vtk, grid, array_name + "_XYZ", vectors, target, names)
        _report(report, vector_name, "DERIVED", target, shape_of(vectors, 3))
        magnitudes = tuple((sqrt(sum(value * value for value in row[:3])),) for row in rows)
        mag_name = add_array(vtk, grid, array_name + "_MAG", magnitudes, target, names)
        _report(report, mag_name, "DERIVED", target, shape_of(magnitudes, 1))
    if "STRESS" in upper and components >= 6:
        rows_mises = tuple((mises(row),) for row in rows)
        legacy_name = add_array(vtk, grid, _mises_name(array_name), rows_mises, target, names)
        _report(report, legacy_name, "DERIVED", target, shape_of(rows_mises, 1))
        if not array_name.endswith("_STRESS"):
            specific_name = add_array(vtk, grid, array_name + "_MISES", rows_mises, target, names)
            _report(report, specific_name, "DERIVED", target, shape_of(rows_mises, 1))


def _indexed_rows(field: Field, count: int, id_to_index: dict[int, int]) -> tuple[tuple[float, ...], ...]:
    rows = [[0.0] * field.cols for _ in range(count)]
    for key, values in field.values.items():
        if isinstance(key, int) and key in id_to_index:
            rows[id_to_index[key]] = list(row_values(values, field.cols))
    return tuple(tuple(row) for row in rows)


def _element_nodal_average(field: Field, ctx: MeshContext) -> tuple[tuple[float, ...], ...]:
    sums = [[0.0] * field.cols for _ in range(ctx.node_count)]
    hits = [0] * ctx.node_count
    for key, values in field.values.items():
        if not isinstance(key, tuple) or len(key) < 2:
            continue
        point_id = ctx.node_to_point.get(key[1])
        if point_id is None:
            continue
        row = row_values(values, field.cols)
        for component, value in enumerate(row):
            sums[point_id][component] += value
        hits[point_id] += 1
    return tuple(tuple(value / hits[index] if hits[index] else 0.0 for value in row) for index, row in enumerate(sums))


def _element_nodal_packed(field: Field, ctx: MeshContext) -> tuple[tuple[float, ...], ...]:
    width = max((len(element.nodes) for element in ctx.elements), default=1) * field.cols
    rows = []
    for element in ctx.elements:
        element_id = required_id(element.id, "element")
        values: list[float] = []
        for node in element.nodes:
            values.extend(row_values(field.values.get((element_id, required_id(node.id, "node")), ()), field.cols))
        values.extend([0.0] * (width - len(values)))
        rows.append(tuple(values))
    return tuple(rows)


def _integration_point_packed(field: Field, ctx: MeshContext) -> tuple[tuple[float, ...], ...]:
    by_element = _by_first_index(field)
    width = max((len(rows) for rows in by_element.values()), default=1) * field.cols
    out = []
    for element in ctx.elements:
        values: list[float] = []
        for row in by_element.get(required_id(element.id, "element"), ()):
            values.extend(row)
        values.extend([0.0] * (width - len(values)))
        out.append(tuple(values))
    return tuple(out)


def _integration_point_average(field: Field, ctx: MeshContext) -> tuple[tuple[float, ...], ...]:
    by_element = _by_first_index(field)
    out = []
    for element in ctx.elements:
        values = by_element.get(required_id(element.id, "element"), ())
        if not values:
            out.append((0.0,) * field.cols)
        else:
            out.append(tuple(sum(row[index] for row in values) / len(values) for index in range(field.cols)))
    return tuple(out)


def _by_first_index(field: Field) -> dict[int, list[tuple[float, ...]]]:
    by_element: dict[int, list[tuple[float, ...]]] = {}
    for key, values in field.values.items():
        if isinstance(key, int):
            element_id = key
        elif isinstance(key, tuple) and key:
            element_id = key[0]
        else:
            continue
        by_element.setdefault(element_id, []).append(row_values(values, field.cols))
    return by_element


def _route_unknown(field: Field, ctx: MeshContext) -> tuple[VtkTarget, tuple[tuple[float, ...], ...]] | None:
    integer_keys = {key for key in field.values if isinstance(key, int)}
    if not integer_keys:
        rows = tuple(row_values(values, field.cols) for _, values in sorted(field.values.items(), key=lambda item: str(item[0])))
        return ("FIELD_DATA", rows) if rows else None
    if integer_keys <= set(ctx.node_to_point):
        return "POINT_DATA", _indexed_rows(field, ctx.node_count, ctx.node_to_point)
    if integer_keys <= set(ctx.element_to_cell):
        return "CELL_DATA", _indexed_rows(field, ctx.cell_count, ctx.element_to_cell)
    return None


def _mises_name(array_name: str) -> str:
    if array_name == "STRESS":
        return "MISES"
    if array_name.endswith("_STRESS"):
        return array_name.removesuffix("_STRESS") + "_MISES"
    return array_name + "_MISES"


def _report(report: ExportReport | None, name: str, action: str, target: str, shape: tuple[int, ...] | None, note: str = "") -> None:
    if report is not None:
        report.add(name, action, target, shape, note)
