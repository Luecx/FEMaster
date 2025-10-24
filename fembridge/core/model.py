from __future__ import annotations

from nodes.nodes import Nodes
from elements.elements import Elements
from constraints.constraints import Constraints

# runtime-light type hints (no hard import dependency at runtime)
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from nodes.node import Node
    from elements.element_base import Element
    from sets.nodeset import NodeSet
    from sets.elementset import ElementSet
    from sets.surfaceset import SurfaceSet         # <-- NEW
    from materials.material_base import Material
    from sections.section_base import Section
    from supports.support_collector import SupportCollector   # <-- NEW
    from loads.load_collector import LoadCollector            # <-- NEW
    from coordinates.coordinate_system_base import CoordinateSystem
    from constraints.constraint_base import ConstraintBase

from sets.nodesets import NodeSets
from sets.elementsets import ElementSets
from sets.surfacesets import SurfaceSets           # <-- NEW
from materials.materials import Materials
from sections.sections import Sections
from profiles.profiles import Profiles
from coordinates.coordinate_systems import CoordinateSystems  # <-- NEW


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

        return "\n".join(parts)

    def to_asami(self) -> str:
        raise NotImplementedError("Model.to_asami not implemented yet.")
