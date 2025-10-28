from __future__ import annotations

import copy

from ..nodes.nodes import Nodes
from ..elements.elements import Elements
from ..constraints.constraints import Constraints
from ..steps.steps import Steps

# runtime-light type hints (no hard import dependency at runtime)
from typing import TYPE_CHECKING, Iterable, Sequence
if TYPE_CHECKING:
    from ..nodes.node import Node
    from ..elements.element_base import Element
    from ..sets.nodeset import NodeSet
    from ..sets.elementset import ElementSet
    from ..sets.surfaceset import SurfaceSet         # <-- NEW
    from ..geometry.segment_group import SegmentGroup
    from ..materials.material_base import Material
    from ..sections.section_base import Section
    from ..supports.support_collector import SupportCollector   # <-- NEW
    from ..loads.load_collector import LoadCollector            # <-- NEW
    from ..coordinates.coordinate_system_base import CoordinateSystem
    from ..constraints.constraint_base import ConstraintBase
    from ..steps.step_base import StepBase

from ..sets.nodesets import NodeSets
from ..sets.elementsets import ElementSets
from ..sets.surfacesets import SurfaceSets           # <-- NEW
from ..materials.materials import Materials
from ..sections.sections import Sections
from ..profiles.profiles import Profiles
from ..coordinates.coordinate_systems import CoordinateSystems  # <-- NEW


class Model:
    """
    Top-level geometry + element + sets + materials + sections container.
    Nodes/Elements/Sets/Sections manage their own id assignment.
    """

    def __init__(self, name: str = "MODEL"):
        self.name = name
        self.nodes = Nodes()
        self.elements = Elements()
        self.node_sets = NodeSets()
        self.element_sets = ElementSets()
        self.surface_sets = SurfaceSets()          # <-- NEW
        self.coordinate_systems = CoordinateSystems()  # <-- NEW
        self.materials = Materials()
        self.profiles = Profiles()
        self.sections = Sections()
        self.elements.model = self
        self.constraints = Constraints()
        self.steps = Steps()

        # --- collectors (simple lists) ---
        self.support_collectors: list["SupportCollector"] = []  # <-- NEW
        self.load_collectors: list["LoadCollector"] = []        # <-- NEW

    # --- adders (delegating to containers) ---
    def add_node(self, node: "Node") -> "Node":
        return self.nodes.add(node)

    def add_element(self, elem: "Element") -> "Element":
        return self.elements.add(elem)

    def add_nodeset(self, s: "NodeSet") -> "NodeSet":
        return self.node_sets.add(s)

    def add_elementset(self, s: "ElementSet") -> "ElementSet":
        return self.element_sets.add(s)

    def add_surfaceset(self, s: "SurfaceSet") -> "SurfaceSet":   # <-- NEW
        return self.surface_sets.add(s)

    def add_material(self, m: "Material") -> "Material":
        return self.materials.add(m)

    def add_section(self, s: "Section") -> "Section":
        return self.sections.add(s)

    def add_profile(self, p):  # <-- existing
        return self.profiles.add(p)

    def add_coordinate_system(self, c: "CoordinateSystem") -> "CoordinateSystem":  # <-- NEW
        return self.coordinate_systems.add(c)

    def add_constraint(self, c: "ConstraintBase") -> "ConstraintBase":
        return self.constraints.add(c)

    # --- collectors adders ---------------------------------------------------
    def add_supportcollector(self, c: "SupportCollector") -> "SupportCollector":  # <-- NEW
        self.support_collectors.append(c)
        return c

    def add_loadcollector(self, c: "LoadCollector") -> "LoadCollector":          # <-- NEW
        self.load_collectors.append(c)
        return c

    # --- step management -------------------------------------------------
    def add_step(self, step: "StepBase") -> "StepBase":
        return self.steps.add(step)

    def extend_steps(self, steps: Iterable["StepBase"]) -> None:
        self.steps.extend(steps)

    def subdivided(self, n: int = 1, *, only_quads: bool = False) -> "Model":
        if n < 1:
            return copy.deepcopy(self)

        element_order = self.elements.determine_element_order()
        if element_order == 0:
            raise ValueError("Element order could not be determined or no elements are present.")
        if element_order == 2:
            raise ValueError("Subdivision is not supported for second-order elements.")

        result = copy.deepcopy(self)
        result.name = f"{self.name}_SUBDIVIDED"
        result.elements.model = result

        for _ in range(n):
            adapter = _SubdivisionAdapter(result)
            edge_midpoints = adapter.compute_edge_midpoints()
            original_elements = list(result.elements._items)
            for element in original_elements:
                if element is None:
                    continue
                adapter.begin_element(element)
                element.subdivide(edge_midpoints, adapter, only_quads=only_quads)

        return result

    def extruded(self, n: int, spacing: float = 1.0) -> "Model":
        from ..nodes.node import Node
        from ..sets.nodeset import NodeSet
        from ..sets.elementset import ElementSet

        if n < 1:
            raise ValueError("Number of extrusion layers n must be >= 1.")

        spacing = float(spacing)
        element_order = self.elements.determine_element_order()
        if element_order == 0:
            raise ValueError("Element order could not be determined or no elements are present.")

        base_node_capacity = len(self.nodes)
        if base_node_capacity == 0:
            raise ValueError("Model has no nodes to extrude.")

        extruded_model = self.__class__(f"{self.name}_EXTRUDED")

        node_lookup: dict[int, Node] = {}
        total_layers = n * element_order + 1
        layer_spacing = spacing / (2.0 if element_order == 2 else 1.0)

        for layer in range(total_layers):
            dz = layer * layer_spacing
            for node in self.nodes._items:
                if node is None or node.node_id is None:
                    continue
                new_node_id = layer * base_node_capacity + node.node_id
                new_node = Node(new_node_id, node.x, node.y, node.z + dz)
                extruded_model.add_node(new_node)
                node_lookup[new_node_id] = new_node

        element_copies: dict[int, list["Element"]] = {}
        for i in range(n):
            for element in self.elements._items:
                if element is None or element.element_id is None:
                    continue

                if element_order == 1:
                    lower = [i * base_node_capacity + nid for nid in element.node_ids]
                    upper = [(i + 1) * base_node_capacity + nid for nid in element.node_ids]
                    node_ids = lower + upper
                else:
                    lower = [i * base_node_capacity + nid for nid in element.node_ids]
                    upper = [(i + 2) * base_node_capacity + nid for nid in element.node_ids]
                    middle = [(i + 1) * base_node_capacity + nid for nid in element.node_ids]
                    half = len(lower) // 2
                    node_ids = (
                        lower[:half]
                        + upper[:half]
                        + lower[half:]
                        + upper[half:]
                        + middle[:half]
                    )

                new_type = f"C3D{len(node_ids)}"
                element_cls = Elements.element_class_for(new_type)
                new_element = extruded_model.add_element(element_cls(None, node_ids))
                element_copies.setdefault(element.element_id, []).append(new_element)

        for node_set in self.node_sets._items:
            if node_set is None:
                continue
            if node_set.name and node_set.name.upper() == "NALL":
                continue

            replicated: list[Node] = []
            for node in node_set:
                nid = getattr(node, "node_id", None)
                if nid is None:
                    continue
                for i in range(n + 1):
                    new_node_id = i * base_node_capacity + nid
                    mapped = node_lookup.get(new_node_id)
                    if mapped is not None:
                        replicated.append(mapped)

            extruded_model.add_nodeset(NodeSet(node_set.name, replicated))

        for elem_set in self.element_sets._items:
            if elem_set is None:
                continue
            if elem_set.name and elem_set.name.upper() == "EALL":
                continue

            replicated: list["Element"] = []
            for elem in elem_set:
                eid = getattr(elem, "element_id", None)
                if eid is None:
                    continue
                replicated.extend(element_copies.get(eid, []))

            extruded_model.add_elementset(ElementSet(elem_set.name, replicated))

        return extruded_model

    @staticmethod
    def mesh_2d(
        segment_groups: Sequence["SegmentGroup"],
        *,
        second_order: bool = False,
        mesh_type: int = 0,
        tolerance: float = 1e-6,
        name: str = "SegmentGroupMesh",
    ) -> "Model":
        import gmsh
        import numpy as np

        from ..nodes.node import Node
        from ..sets.nodeset import NodeSet

        if not segment_groups:
            raise ValueError("mesh_2d requires at least one SegmentGroup.")

        gmsh.initialize()
        try:
            gmsh.model.add(name or "SegmentGroupMesh")
            gmsh.option.setNumber("General.Terminal", 0)

            curve_loops: list[int] = []
            internal_curves: list[int] = []

            for group in segment_groups:
                points = group.get_points()
                is_closed = group.is_closed()
                if is_closed:
                    points = points[:-1]

                point_tags = [gmsh.model.occ.addPoint(float(x), float(y), 0.0) for x, y in points]

                if is_closed:
                    n = len(point_tags)
                    line_tags = [gmsh.model.occ.addLine(point_tags[i], point_tags[(i + 1) % n]) for i in range(n)]
                    curve_loops.append(gmsh.model.occ.addCurveLoop(line_tags))
                else:
                    for a, b in zip(point_tags[:-1], point_tags[1:]):
                        internal_curves.append(gmsh.model.occ.addLine(a, b))

            if not curve_loops:
                raise ValueError("At least one closed SegmentGroup is required to create a meshable surface.")

            surface_tag = gmsh.model.occ.addPlaneSurface(curve_loops)
            gmsh.model.occ.synchronize()

            if internal_curves:
                gmsh.model.occ.fragment(
                    [(2, surface_tag)],
                    [(1, tag) for tag in internal_curves],
                )
                gmsh.model.occ.synchronize()

            gmsh.option.setNumber("Mesh.Algorithm", 6)
            gmsh.option.setNumber("Mesh.ElementOrder", 2 if second_order else 1)
            if mesh_type in {1, 2}:
                gmsh.option.setNumber("Mesh.RecombineAll", 1)

            gmsh.model.mesh.generate(2)

            node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
            coords = node_coords.reshape((-1, 3))

            mesh_model = Model(name)
            tag_to_node: dict[int, "Node"] = {}
            for tag, (x, y, _) in zip(node_tags, coords):
                node = mesh_model.add_node(Node(None, x, y, 0.0))
                tag_to_node[int(tag)] = node

            elt_types, _, elt_node_tags = gmsh.model.mesh.getElements(dim=2)
            node_counts = {2: 3, 3: 4, 9: 6, 10: 9}
            type_lookup = {3: "S3", 4: "S4", 6: "S6", 8: "S8"}

            for etype, flat_nodes in zip(elt_types, elt_node_tags):
                count = node_counts.get(int(etype))
                if count is None:
                    continue

                nodes_array = np.array(flat_nodes, dtype=int).reshape(-1, count)
                for row in nodes_array:
                    node_objs = [tag_to_node[int(tag)] for tag in row]
                    if count == 9:
                        node_objs = node_objs[:8]
                    elem_type = type_lookup[len(node_objs)]
                    node_ids = [node.node_id for node in node_objs if node.node_id is not None]
                    element_cls = Elements.element_class_for(elem_type)
                    mesh_model.add_element(element_cls(None, node_ids))

            group_sets: dict[int, NodeSet] = {}
            segment_sets: dict[int, NodeSet] = {}
            for group in segment_groups:
                group_set = NodeSet(group.name, [])
                mesh_model.add_nodeset(group_set)
                group_sets[id(group)] = group_set
                for segment in group.segments:
                    segment_set = NodeSet(segment.name, [])
                    mesh_model.add_nodeset(segment_set)
                    segment_sets[id(segment)] = segment_set

            for node in mesh_model.nodes._items:
                if node is None:
                    continue
                point = np.array([node.x, node.y])
                for group in segment_groups:
                    group_set = group_sets[id(group)]
                    for segment in group.segments:
                        if segment.contains(point, tolerance=tolerance):
                            segment_sets[id(segment)].add(node)
                            group_set.add(node)

            return mesh_model
        finally:
            if gmsh.isInitialized():
                gmsh.finalize()

    # --- serialization ---
    def to_femaster(self) -> str:
        parts = [
            f"*MODEL, NAME={self.name}",
            self.nodes.to_femaster(),
            self.elements.to_femaster(),
            self.node_sets.to_femaster(),
            self.element_sets.to_femaster(),
            self.surface_sets.to_femaster(),
            self.materials.to_femaster(),
            self.profiles.to_femaster(),
            self.sections.to_femaster(),
            self.coordinate_systems.to_femaster(),
        ]

        constraints_block = self.constraints.to_femaster()
        if constraints_block:
            parts.append(constraints_block)

        # append collectors: supports first, then loads
        for sc in self.support_collectors:
            parts.append(sc.to_femaster())
        for lc in self.load_collectors:
            parts.append(lc.to_femaster())

        steps_block = self.steps.to_femaster()
        if steps_block:
            parts.append(steps_block)

        return "\n".join(parts)

    def to_asami(self) -> str:
        raise NotImplementedError("Model.to_asami not implemented yet.")


class _SubdivisionAdapter:
    def __init__(self, model: "Model") -> None:
        self.model = model
        self.nodes = model.nodes
        self.elements = model.elements
        self.membership: dict[int, list[tuple["ElementSet", int]]] = {}
        self._build_membership()
        self.current_parent_id: int | None = None
        self.current_positions: list[tuple["ElementSet", int]] = []

    def _build_membership(self) -> None:
        from ..sets.elementset import ElementSet

        for elem_set in self.model.element_sets._items:
            if elem_set is None:
                continue
            if not isinstance(elem_set, ElementSet):
                continue
            for idx, element in enumerate(elem_set.elements):
                eid = getattr(element, "element_id", None)
                if eid is None:
                    continue
                self.membership.setdefault(eid, []).append((elem_set, idx))

    def begin_element(self, element: "Element") -> None:
        eid = getattr(element, "element_id", None)
        if eid is None:
            self.current_parent_id = None
            self.current_positions = []
            return
        self.current_parent_id = eid
        self.current_positions = list(self.membership.get(eid, []))

    def add_node(self, node_id: int | None = None, *, x: float, y: float, z: float) -> int:
        from ..nodes.node import Node

        new_node = Node(node_id, x, y, z)
        added = self.model.add_node(new_node)
        return int(added.node_id)

    def add_element(
        self,
        *,
        element_type: str,
        node_ids: Iterable[int],
        element_id: int | None = None,
    ) -> int:
        node_list = list(node_ids)
        element_cls = Elements.element_class_for(element_type)
        new_element = self.model.add_element(element_cls(element_id, node_list))
        eid = int(new_element.element_id)

        if element_id is not None:
            positions = self.membership.get(element_id, [])
            if positions:
                for elem_set, idx in positions:
                    elem_set.elements[idx] = new_element
                self.membership[element_id] = [(elem_set, idx) for elem_set, idx in positions]
            else:
                self.membership[element_id] = []
        else:
            entries: list[tuple["ElementSet", int]] = []
            for elem_set, _ in self.current_positions:
                elem_set.add(new_element)
                entries.append((elem_set, len(elem_set.elements) - 1))
            if entries:
                self.membership[eid] = entries
        return eid

    def compute_edge_midpoints(self) -> dict[tuple[int, int], int]:
        edge_midpoints: dict[tuple[int, int], int] = {}
        for element in self.elements._items:
            if element is None:
                continue
            for n1, n2 in element.connectivity():
                if (n1, n2) in edge_midpoints:
                    continue
                node1 = self.nodes[n1]
                node2 = self.nodes[n2]
                x = 0.5 * (node1.x + node2.x)
                y = 0.5 * (node1.y + node2.y)
                z = 0.5 * (node1.z + node2.z)
                new_id = self.add_node(x=x, y=y, z=z)
                edge_midpoints[(n1, n2)] = new_id
                edge_midpoints[(n2, n1)] = new_id
        return edge_midpoints
