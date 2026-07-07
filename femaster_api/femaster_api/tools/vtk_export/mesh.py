"""Mesh, set, section, and coupling export helpers."""

from __future__ import annotations

from femaster_api.model import EntityType, Model
from femaster_api.model.constraints import CouplingConstraint, CouplingType
from femaster_api.model.sections import ShellSection

from .common import MeshContext, NameRegistry, add_array, add_string_array, cell_type, required_id


def build_grid(vtk, model: Model, *, export_sets: bool = False, pack_sets: bool = False, export_shell_thickness: bool = False, export_couplings: bool = False):
    nodes = tuple(model.nodes)
    elements = tuple(model.elements)
    node_to_point = {required_id(node.id, "node"): index for index, node in enumerate(nodes)}
    element_to_cell: dict[int, int] = {}

    grid = vtk.vtkUnstructuredGrid()
    points = vtk.vtkPoints()
    points.SetNumberOfPoints(len(nodes))
    for index, node in enumerate(nodes):
        points.SetPoint(index, float(node.x), float(node.y), float(node.z))
    grid.SetPoints(points)

    for element in elements:
        ids = vtk.vtkIdList()
        for node in element.nodes:
            ids.InsertNextId(node_to_point[required_id(node.id, "node")])
        cell_id = grid.GetNumberOfCells()
        grid.InsertNextCell(cell_type(vtk, element), ids)
        element_to_cell[required_id(element.id, "element")] = cell_id

    names = NameRegistry()
    if export_couplings:
        _add_coupling_cells(vtk, grid, model, node_to_point, names)

    ctx = MeshContext(nodes, elements, node_to_point, element_to_cell, grid.GetNumberOfCells())

    if export_shell_thickness:
        add_array(vtk, grid, "SHELL_THICKNESS", shell_thickness_rows(model, ctx), "CELL_DATA", names)

    if export_sets:
        add_sets(vtk, grid, model, ctx, names, pack=pack_sets)

    return grid, ctx, names


def add_sets(vtk, grid, model: Model, ctx: MeshContext, names: NameRegistry, *, pack: bool) -> None:
    node_sets = tuple(item for item in model.sets.all(EntityType.NODE) if item.name)
    element_sets = tuple(item for item in model.sets.all(EntityType.ELEMENT) if item.name)
    if pack:
        if node_sets:
            add_array(vtk, grid, "SETS_N", _packed_node_sets(node_sets, ctx), "POINT_DATA", names, integer=True)
            add_string_array(vtk, grid, "SETS_N_NAMES", [item.name for item in node_sets], names)
        if element_sets:
            add_array(vtk, grid, "SETS_E", _packed_element_sets(element_sets, ctx), "CELL_DATA", names, integer=True)
            add_string_array(vtk, grid, "SETS_E_NAMES", [item.name for item in element_sets], names)
        return
    for item in node_sets:
        add_array(vtk, grid, f"SET_N_{item.name}", _single_node_set(item, ctx), "POINT_DATA", names, integer=True)
    for item in element_sets:
        add_array(vtk, grid, f"SET_E_{item.name}", _single_element_set(item, ctx), "CELL_DATA", names, integer=True)


def shell_thickness_rows(model: Model, ctx: MeshContext) -> tuple[tuple[float, ...], ...]:
    rows = [[0.0] for _ in range(ctx.cell_count)]
    for section in model.sections:
        if not isinstance(section, ShellSection):
            continue
        for element in section.element_set.members:
            cell_id = ctx.element_to_cell.get(required_id(element.id, "element"))
            if cell_id is not None:
                rows[cell_id][0] = float(section.thickness)
    return tuple(tuple(row) for row in rows)


def _add_coupling_cells(vtk, grid, model: Model, node_to_point: dict[int, int], names: NameRegistry) -> None:
    role = [[0] for _ in range(grid.GetNumberOfCells())]
    coupling_id = [[0] for _ in range(grid.GetNumberOfCells())]
    dofs = [[0, 0, 0, 0, 0, 0] for _ in range(grid.GetNumberOfCells())]
    counter = 0

    for constraint in model.constraints:
        if not isinstance(constraint, CouplingConstraint) or constraint.type is not CouplingType.KINEMATIC:
            continue
        master_nodes = _nodes_from_region(constraint.master)
        slave_nodes = _nodes_from_region(constraint.slave)
        if not master_nodes or not slave_nodes:
            continue
        counter += 1
        master_pid = node_to_point.get(required_id(master_nodes[0].id, "node"))
        if master_pid is None:
            continue
        dof_row = [int(value) for value in constraint.dofs]
        for slave in slave_nodes:
            slave_pid = node_to_point.get(required_id(slave.id, "node"))
            if slave_pid is None:
                continue
            line = vtk.vtkLine()
            ids = line.GetPointIds()
            ids.SetId(0, master_pid)
            ids.SetId(1, slave_pid)
            grid.InsertNextCell(line.GetCellType(), ids)
            role.append([1])
            coupling_id.append([counter])
            dofs.append(dof_row)

    if counter:
        add_array(vtk, grid, "CELL_ROLE", tuple(tuple(row) for row in role), "CELL_DATA", names, integer=True)
        add_array(vtk, grid, "COUPLING_ID", tuple(tuple(row) for row in coupling_id), "CELL_DATA", names, integer=True)
        add_array(vtk, grid, "COUPLING_DOF_MASK", tuple(tuple(row) for row in dofs), "CELL_DATA", names, integer=True)


def _nodes_from_region(region) -> tuple[object, ...]:
    entity_type = getattr(region, "entity_type", None)
    if entity_type == EntityType.NODE:
        return tuple(region.members)
    if entity_type == EntityType.ELEMENT:
        seen: dict[int, object] = {}
        for element in region.members:
            for node in getattr(element, "nodes", ()):
                if node.id is not None:
                    seen[int(node.id)] = node
        return tuple(seen[key] for key in sorted(seen))
    return ()


def _single_node_set(node_set, ctx: MeshContext) -> tuple[tuple[float, ...], ...]:
    rows = [[0] for _ in range(ctx.node_count)]
    for node in node_set.members:
        point_id = ctx.node_to_point.get(required_id(node.id, "node"))
        if point_id is not None:
            rows[point_id][0] = 1
    return tuple(tuple(row) for row in rows)


def _single_element_set(element_set, ctx: MeshContext) -> tuple[tuple[float, ...], ...]:
    rows = [[0] for _ in range(ctx.cell_count)]
    for element in element_set.members:
        cell_id = ctx.element_to_cell.get(required_id(element.id, "element"))
        if cell_id is not None:
            rows[cell_id][0] = 1
    return tuple(tuple(row) for row in rows)


def _packed_node_sets(node_sets, ctx: MeshContext) -> tuple[tuple[float, ...], ...]:
    rows = [[0] * len(node_sets) for _ in range(ctx.node_count)]
    for set_id, node_set in enumerate(node_sets):
        for node in node_set.members:
            point_id = ctx.node_to_point.get(required_id(node.id, "node"))
            if point_id is not None:
                rows[point_id][set_id] = 1
    return tuple(tuple(row) for row in rows)


def _packed_element_sets(element_sets, ctx: MeshContext) -> tuple[tuple[float, ...], ...]:
    rows = [[0] * len(element_sets) for _ in range(ctx.cell_count)]
    for set_id, element_set in enumerate(element_sets):
        for element in element_set.members:
            cell_id = ctx.element_to_cell.get(required_id(element.id, "element"))
            if cell_id is not None:
                rows[cell_id][set_id] = 1
    return tuple(tuple(row) for row in rows)
