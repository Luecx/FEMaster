from .fem_elements import *

from .fem_geometry_mesh_2d import mesh_interior
from .fem_geometry_subdivided import subdivided_geometry
from .fem_geometry_extruded import extruded_geometry
from .fem_geometry_write_inp import write_input_deck
from .fem_geometry_read_inp import read_input_deck
from .fem_geometry_connectivity import connectivity_node_to_element
from .fem_geometry_connectivity import connectivity_element_to_element
from .fem_geometry_connectivity import element_element_distance_matrix

import re
import numpy as np
import copy
import scipy.sparse
import matplotlib.pyplot as plt

from collections import Counter


class Geometry:

    element_classes = {
        'C2D3': C2D3,
        'C2D4': C2D4,
        'C2D6': C2D6,
        'C2D8': C2D8,
        'C3D4': C3D4,
        'C3D5': C3D5,
        'C3D6': C3D6,
        'C3D8': C3D8,
        'C3D10': C3D10,
        'C3D15': C3D15,
        'C3D20': C3D20,
        'C3D20R': C3D20,
        'B33': B33,
        'S4': S4,
        'S8': S8,
        'S3': S3,
        'S6': S6,
        'P': Point,
    }

    def __init__(self, dimension=3):
        self.dimension = dimension  # 2 for 2D, 3 for 3D
        self.nodes = []
        self.elements = []
        self.node_sets = {"NALL": []}
        self.elem_sets = {"EALL": []}

        # Materialien / Sektionen / Profile
        self.materials = {}
        self.sections  = []
        self.profiles  = []

        # Collector (Name -> {set, data, coord_sys})
        # data = (fx, fy, fz, mx, my, mz)
        self.loads     = {}
        self.supps     = {}

        # Steps (Liste): z.B. {'name':..., 'type':'LINEAR STATIC', 'loads':[...], 'supps':[...]}
        self.steps     = []

        # Koordinatensysteme (nur Rohdaten speichern, wie im Input angegeben)
        # name -> {'type': 'RECTANGULAR'|'CYLINDRICAL',
        #          'values': tuple[float, ...]}  # 3/6/9 bei RECT, 9 bei CYL
        # Optional kann zusätzlich 'definition' für Schreibzwecke stehen (z.B. 'VECTOR' oder 'POINTS').
        self.coordinate_systems = {}

        # Couplings (Liste reiner Daten fürs Input-Deck)
        # {'type': 'KINEMATIC', 'master': 'M_*', 'slave': 'C_*',
        #  'cx':int, 'cy':int, 'cz':int, 'crx':int, 'cry':int, 'crz':int}
        self.couplings = []

        # Connector-Elements (Liste reiner Daten fürs Input-Deck)
        # {'type':'HINGE'|'BEAM'|'CYLINDRICAL'|'TRANSLATOR'|'JOIN',
        #  'coord_sys':'CSY', 'nset1':'M1', 'nset2':'M2'}
        self.connectors = []

    # -------------------------
    # Nodes / Elements / Sets
    # -------------------------
    def add_node(self, node_id=-1, x=0, y=0, z=0):
        # Resize nodes list if necessary
        if node_id == -1:
            node_id = len(self.nodes)
        if node_id >= len(self.nodes):
            self.nodes.extend([None] * (node_id - len(self.nodes) + 1))
        self.nodes[node_id] = [x, y, z]
        self.node_sets["NALL"].append(node_id)
        return node_id

    def add_element(self, element_id=-1, element_type=None, node_ids=None):
        # Create the appropriate element type
        if element_type not in Geometry.element_classes:
            raise ValueError(f"Unsupported element type: {element_type}")

        if element_id == -1:
            element_id = len(self.elements)

        element = Geometry.element_classes[element_type](element_id, node_ids)

        # Resize elements list if necessary
        if element_id >= len(self.elements):
            self.elements.extend([None] * (element_id - len(self.elements) + 1))
        self.elements[element_id] = element
        self.elem_sets["EALL"].append(element_id)
        return element_id

    def add_element_beams(self,
                          node1: int,
                          node2: int,
                          element_id: int | None = None,
                          orientation: tuple[float, float, float] | list[float] | None = None,
                          *,
                          material: str | None = None,
                          profile: str | None = None,
                          n_elements: int = 1,
                          joint_A: bool = False,
                          joint_B: bool = False):
        """
        Create one or several B33 beam elements between node1 (A) and node2 (B).
        Returns int if n_elements == 1 else list[int].
        Uses NumPy for all vector math.
        """
        # --- basic check
        if node1 >= len(self.nodes) or self.nodes[node1] is None:
            raise ValueError(f"node1 id {node1} does not exist")
        if node2 >= len(self.nodes) or self.nodes[node2] is None:
            raise ValueError(f"node2 id {node2} does not exist")
        if n_elements < 1:
            raise ValueError("n_elements must be >= 1")

        # --- ensure a rectangular CSYS for connectors
        def _ensure_default_rect_csys() -> str:
            for name, cs in self.coordinate_systems.items():
                if cs.get('type') == 'RECTANGULAR':
                    return name
            name = "CSYS_GLOBAL"
            self.add_coordinate_system_rectangular(name, 1.0, 0.0, 0.0, definition="VECTOR")
            return name

        csys_name = _ensure_default_rect_csys()

        # --- joints (only on requested sides)
        A_id = node1
        B_id = node2

        def _make_joint(side_label: str, old_id: int) -> int:
            xyz    = self.nodes[old_id]
            new_id = self.add_node(x=xyz[0], y=xyz[1], z=xyz[2])
            set_new = f"JOINT_{side_label}_NEW_{new_id}"
            set_old = f"JOINT_{side_label}_OLD_{new_id}"
            self.add_node_set(set_new)
            self.add_node_to_set(set_new, new_id)
            self.add_node_set(set_old)
            self.add_node_to_set(set_old, old_id)
            self.add_connector(type="JOIN", coord_sys=csys_name, nset1=set_new, nset2=set_old)
            return new_id

        if joint_A:
            A_id = _make_joint("A", A_id)
        if joint_B:
            B_id = _make_joint("B", B_id)

        # --- NumPy helpers
        def _np_norm(v: np.ndarray) -> float:
            return float(np.linalg.norm(v))

        def _normalize(v: np.ndarray) -> np.ndarray:
            n = _np_norm(v)
            return v / n if n > 1e-14 else np.zeros(3)

        def _auto_orientation(a_xyz: np.ndarray, b_xyz: np.ndarray) -> np.ndarray:
            d = _normalize(b_xyz - a_xyz)                       # beam axis
            ez = np.array([0.0, 0.0, 1.0])
            ex = np.array([1.0, 0.0, 0.0])
            helper = ez if abs(float(np.dot(d, ez))) < 0.9 else ex
            # Gram-Schmidt: helper ⟂ d
            v = helper - np.dot(helper, d) * d
            v = _normalize(v)
            if _np_norm(v) <= 1e-14:
                # fallback: any perpendicular via cross
                v = _normalize(np.cross(d, ex))
            return v

        A_xyz = np.asarray(self.nodes[A_id], dtype=float)
        B_xyz = np.asarray(self.nodes[B_id], dtype=float)

        # --- orientation: param or auto
        if orientation is not None:
            orientation = _normalize(np.asarray(orientation, dtype=float))
            if _np_norm(orientation) <= 1e-14:
                orientation = _auto_orientation(A_xyz, B_xyz)
        else:
            orientation = _auto_orientation(A_xyz, B_xyz)

        # --- choose material/profile (auto if missing)
        if material is None:
            if not self.materials:
                raise ValueError("No materials defined; please add a material or pass 'material'.")
            material = next(iter(self.materials.keys()))
        if profile is None:
            if not self.profiles:
                raise ValueError("No profiles defined; please add a profile or pass 'profile'.")
            first_prof = self.profiles[0]
            profile = first_prof['name'] if isinstance(first_prof, dict) and 'name' in first_prof else first_prof

        # --- subdivision nodes with NumPy (linspace)
        path_nodes = [A_id]
        if n_elements > 1:
            # internal points at k/n_elements, k = 1..n_elements-1
            ts = np.linspace(0.0, 1.0, n_elements + 1)[1:-1]
            seg = (B_xyz - A_xyz)
            for t in ts:
                p = A_xyz + t * seg
                nid = self.add_node(x=float(p[0]), y=float(p[1]), z=float(p[2]))
                path_nodes.append(nid)
        path_nodes.append(B_id)

        # --- element id policy
        fixed_ids = None
        if element_id is not None:
            fixed_ids = [element_id + i for i in range(n_elements)]
            for eid in fixed_ids:
                if eid < len(self.elements) and self.elements[eid] is not None:
                    raise ValueError(f"Requested element_id {eid} is already occupied.")

        # --- create elements
        created = []
        for i in range(n_elements):
            nA = path_nodes[i]
            nB = path_nodes[i + 1]
            if fixed_ids is not None:
                eid = self.add_element(element_id=fixed_ids[i], element_type='B33', node_ids=[nA, nB])
            else:
                eid = self.add_element(element_id=-1, element_type='B33', node_ids=[nA, nB])
            created.append(eid)

        # --- attach a BEAM section for these elements (new ELSET)
        elset_name = f"EL_BEAM_{created[0]}_{created[-1]}"
        self.add_element_set(elset_name)
        for eid in created:
            self.add_element_to_set(elset_name, eid)

        self.add_section_beam(elset=elset_name,
                              material=material,
                              profile=profile,
                              orientation=(float(orientation[0]),
                                           float(orientation[1]),
                                           float(orientation[2])))

        return created[0] if n_elements == 1 else created


    def add_node_set(self, name):
        if name not in self.node_sets:
            self.node_sets[name] = []

    def add_node_to_set(self, name, node_id):
        self.add_node_set(name)
        self.node_sets[name].append(node_id)

    def add_element_set(self, name):
        if name not in self.elem_sets:
            self.elem_sets[name] = []

    def add_element_to_set(self, name, element_id):
        self.add_element_set(name)
        self.elem_sets[name].append(element_id)

    # -------------------------
    # Materials / Sections / Profiles
    # -------------------------
    def add_material(self, name, young, poisson, density):
        self.materials[name] = {'young': young, 'poisson': poisson, 'density': density}

    def add_profile(self, name, area, Iy, Iz, It):
        self.profiles.append({
            'name': name,
            'area': area,
            'Iy': Iy,
            'Iz': Iz,
            'It': It,
        })

    def add_section_solid(self, material, elset):
        self.sections.append({'type': 'SOLID', 'material': material, 'elset': elset})
        return self.sections[-1]

    def add_section_shell(self, material, elset, thickness):
        self.sections.append({'type': 'SHELL', 'material': material, 'elset': elset, 'thickness': thickness})
        return self.sections[-1]

    def add_section_beam(self, elset, material, profile, orientation=(0, 1, 0)):
        if not isinstance(orientation, (list, tuple)) or len(orientation) != 3:
            raise ValueError("Orientation must be a list or tuple of 3 floats")
        self.sections.append({
            'type': 'BEAM',
            'elset': elset,
            'material': material,
            'profile': profile,
            'orientation': tuple(orientation)
        })
        return self.sections[-1]

    def add_section_pointmass(self, elset, mass=None, inertia=None, spring=None, rotary_spring=None):
        def validate_vector(v, name):
            if v is not None:
                if not isinstance(v, (list, tuple)) or len(v) != 3:
                    raise ValueError(f"{name} must be a list of 3 values or None.")
                return list(v)
            return None

        inertia = validate_vector(inertia, "inertia")
        spring = validate_vector(spring, "spring")
        rotary_spring = validate_vector(rotary_spring, "rotary_spring")

        self.sections.append({
            'type': 'POINTMASS',
            'elset': elset,
            'mass': mass,
            'inertia': inertia,
            'spring': spring,
            'rotary_spring': rotary_spring
        })

    # -------------------------
    # Coordinate Systems (nur Rohdaten speichern)
    # -------------------------
    def add_coordinate_system_rectangular(self, name, *values, definition="VECTOR"):
        """
        Speichert ein RECTANGULAR-System, genau wie im Input angegeben.
        Erlaubt 3, 6 oder 9 Werte (x, y [, z]).
        values: (x1,y1,z1[, x2,y2,z2[, x3,y3,z3]])
        """
        if name in self.coordinate_systems:
            raise ValueError(f"Coordinate system '{name}' already exists.")
        if len(values) not in (3, 6, 9):
            raise ValueError("RECTANGULAR coordinate system requires 3, 6, or 9 values.")
        self.coordinate_systems[name] = {
            'type': 'RECTANGULAR',
            'definition': definition,
            'values': tuple(float(v) for v in values)
        }
        return name

    def add_coordinate_system_cylindrical(self, name, *values, definition="POINTS"):
        """
        Speichert ein CYLINDRICAL-System, genau wie im Input angegeben.
        Erlaubt 9 Werte (Base, Radial, Theta-Punkt): (x1,y1,z1, x2,y2,z2, x3,y3,z3)
        """
        if name in self.coordinate_systems:
            raise ValueError(f"Coordinate system '{name}' already exists.")
        if len(values) != 9:
            raise ValueError("CYLINDRICAL coordinate system requires 9 values (base, radial, theta points).")
        self.coordinate_systems[name] = {
            'type': 'CYLINDRICAL',
            'definition': definition,
            'values': tuple(float(v) for v in values)
        }
        return name

    # -------------------------
    # Loads / Supports (Collectors) – optionales CSYS
    # -------------------------
    def add_load(self, name, set, fx=None, fy=None, fz=None, mx=None, my=None, mz=None, coord_sys=None):
        if coord_sys is not None and coord_sys not in self.coordinate_systems:
            raise KeyError(f"Unknown coordinate system '{coord_sys}' for load '{name}'.")
        self.loads[name] = {'set': set, 'data': (fx, fy, fz, mx, my, mz), 'coord_sys': coord_sys}

    def add_supp(self, name, set, fx=None, fy=None, fz=None, mx=None, my=None, mz=None, coord_sys=None):
        if coord_sys is not None and coord_sys not in self.coordinate_systems:
            raise KeyError(f"Unknown coordinate system '{coord_sys}' for support '{name}'.")
        self.supps[name] = {'set': set, 'data': (fx, fy, fz, mx, my, mz), 'coord_sys': coord_sys}

    # -------------------------
    # Steps
    # -------------------------
    def add_step(self, type, name="Step-1", **kwargs):
        """
        Aktuelle Step-Typen:
          - "LINEAR STATIC": required keys = ['loads', 'supps']
          - "EIGENFREQ"    : required keys = ['supps', 'numeigenvalues']
        """
        # remove spaces from type
        type = type.replace(" ", "").upper()

        step = {'name': name, 'type': type}

        if step['type'] == "LINEARSTATIC":
            required_keys = ['loads', 'supps']
        elif step['type'] == "EIGENFREQ":
            required_keys = ['supps', 'numeigenvalues']
        elif step['type'] == "LINEARBUCKLING":
            required_keys = ['loads', 'supps', 'numeigenvalues']
        else:
            raise ValueError(f"Unsupported step type: {type}")

        for key in required_keys:
            if key not in kwargs:
                raise ValueError(f"{key} must be provided for step type '{type}'")

        step.update(kwargs)
        self.steps.append(step)

    # -------------------------
    # Couplings (explizite DOFs)
    # -------------------------
    def add_coupling(self, master_set, slave_set, type="KINEMATIC",
                     cx=1, cy=1, cz=1, crx=1, cry=1, crz=1):
        """
        Speichert einen Coupling-Eintrag für das Input-Deck.
        master_set: Name eines Node-Sets mit GENAU 1 Knoten (Referenz).
        slave_set : Name eines Node-Sets (typisch mehrere Knoten).
        type      : aktuell 'KINEMATIC' (weitere können später ergänzt werden).
        DOFs      : cx, cy, cz, crx, cry, crz (0/1).
        """
        t = str(type).upper()
        if t != "KINEMATIC":
            raise ValueError(f"Unsupported coupling type '{type}'. Only KINEMATIC is implemented currently.")

        if master_set not in self.node_sets:
            raise KeyError(f"Unknown node set '{master_set}' (master).")
        if slave_set not in self.node_sets:
            raise KeyError(f"Unknown node set '{slave_set}' (slave).")

        m_nodes = self.node_sets[master_set]
        if len(m_nodes) != 1:
            raise ValueError(f"Master node set '{master_set}' must contain exactly one node (got {len(m_nodes)}).")

        self.couplings.append({
            'type': t,
            'master': master_set,
            'slave': slave_set,
            'cx': int(bool(cx)),
            'cy': int(bool(cy)),
            'cz': int(bool(cz)),
            'crx': int(bool(crx)),
            'cry': int(bool(cry)),
            'crz': int(bool(crz)),
        })

    # -------------------------
    # Connector Elements (nur Rohdaten)
    # -------------------------
    def add_connector(self, type, coord_sys, nset1, nset2):
        """
        type     : 'BEAM'|'HINGE'|'CYLINDRICAL'|'TRANSLATOR'|'JOIN' (aktuell bekannte)
        coord_sys: Name eines existierenden Koordinatensystems (wird nur referenziert)
        nset1/2  : Node-Set-Namen, JEWEILS mit GENAU 1 Knoten.
        """
        t = str(type).upper()
        if t not in ['BEAM', 'HINGE', 'CYLINDRICAL', 'TRANSLATOR', 'JOIN']:
            raise ValueError(f"Unsupported connector type '{type}'.")
        if coord_sys not in self.coordinate_systems:
            raise KeyError(f"Unknown coordinate system '{coord_sys}' for connector.")
        if nset1 not in self.node_sets:
            raise KeyError(f"Unknown node set '{nset1}' for connector.")
        if nset2 not in self.node_sets:
            raise KeyError(f"Unknown node set '{nset2}' for connector.")
        if len(self.node_sets[nset1]) != 1:
            raise ValueError(f"Connector nset1 '{nset1}' must contain exactly one node.")
        if len(self.node_sets[nset2]) != 1:
            raise ValueError(f"Connector nset2 '{nset2}' must contain exactly one node.")

        self.connectors.append({
            'type': t,
            'coord_sys': coord_sys,
            'nset1': nset1,
            'nset2': nset2
        })

    # -------------------------
    # Ableitungen / Utilities
    # -------------------------
    def determine_element_order(self):
        for element in self.elements:
            if element is not None:
                if element.elem_type in ['C2D3', 'C2D4', 'C3D4', 'C3D6', 'C3D8']:
                    return 1
                if element.elem_type in ['S3', 'S4', 'S6', 'S8']:
                    return 1
                elif element.elem_type in ['C2D6', 'C2D8', 'C3D10', 'C3D15', 'C3D20']:
                    return 2
        return 0

    def get_edges(self, sort=True):
        edges = set()
        for element in self.elements:
            if element is not None:
                for edge in element.connectivity():
                    if sort:
                        edges |= {tuple(sorted(edge))}
                    else:
                        edges |= {edge}
        return edges

    def compute_edge_midpoints(self):
        edge_midpoints = {}
        # Collect all edges and create midpoints
        for element in self.elements:
            if element is not None:
                for edge in element.connectivity():

                    id1 = edge[0]
                    id2 = edge[1]

                    if (id1, id2) not in edge_midpoints:
                        n1 = self.nodes[id1]
                        n2 = self.nodes[id2]
                        midpoint = [(n1[i] + n2[i]) / 2 for i in range(3)]
                        new_node_id = len(self.nodes)
                        self.add_node(new_node_id, *midpoint)

                        # store the new node id for both directions
                        edge_midpoints[(id1, id2)] = new_node_id
                        edge_midpoints[(id2, id1)] = new_node_id
        return edge_midpoints

    def compute_element_midpoints(self):
        midpoints = [(0, 0, 0) for _ in range(len(self.elements))]
        for i, element in enumerate(self.elements):
            if element is None:
                continue

            element_node_ids = element.node_ids
            element_nodes = [self.nodes[node_id] for node_id in element_node_ids]
            center = np.mean(element_nodes, axis=0)
            midpoints[i] = center
        return midpoints

    def as_second_order(self):
        # Create a new Geometry object for the second-order elements
        geometry = Geometry(self.dimension)

        # Copy nodes and node sets
        geometry.nodes = copy.deepcopy(self.nodes)
        geometry.node_sets = {name: ids.copy() for name, ids in self.node_sets.items()}
        geometry.elem_sets = {name: ids.copy() for name, ids in self.elem_sets.items()}
        geometry.elements = copy.deepcopy(self.elements)

        # Kopiere sonstige Strukturen 1:1 (ohne Änderungen)
        geometry.materials = copy.deepcopy(self.materials)
        geometry.sections = copy.deepcopy(self.sections)
        geometry.profiles = copy.deepcopy(self.profiles)
        geometry.loads = copy.deepcopy(self.loads)
        geometry.supps = copy.deepcopy(self.supps)
        geometry.steps = copy.deepcopy(self.steps)
        geometry.coordinate_systems = copy.deepcopy(self.coordinate_systems)  # NICHT spiegeln/ändern
        geometry.couplings = copy.deepcopy(self.couplings)
        geometry.connectors = copy.deepcopy(self.connectors)

        edge_midpoints = geometry.compute_edge_midpoints()

        # Convert elements to second order and replace inplace
        for i in range(len(self.elements)):
            element = self.elements[i]
            if element is not None:
                geometry.elements[i] = element.to_second_order(edge_midpoints)
        return geometry

    def mirrored(self, plane, location, copy_nodes=True):
        geometry = Geometry(self.dimension)
        geometry.nodes = copy.deepcopy(self.nodes)
        geometry.node_sets = {name: ids.copy() for name, ids in self.node_sets.items()}
        geometry.elem_sets = {name: ids.copy() for name, ids in self.elem_sets.items()}
        geometry.elements = copy.deepcopy(self.elements)

        # Alles andere nur kopieren, nicht transformieren
        geometry.materials = copy.deepcopy(self.materials)
        geometry.sections = copy.deepcopy(self.sections)
        geometry.profiles = copy.deepcopy(self.profiles)
        geometry.loads = copy.deepcopy(self.loads)
        geometry.supps = copy.deepcopy(self.supps)
        geometry.steps = copy.deepcopy(self.steps)
        geometry.coordinate_systems = copy.deepcopy(self.coordinate_systems)  # unverändert lassen
        geometry.couplings = copy.deepcopy(self.couplings)
        geometry.connectors = copy.deepcopy(self.connectors)

        axis_idx = ['yz', 'xz', 'xy'].index(plane)

        if not copy_nodes:
            # mirror all nodes and elements inplace
            for node in geometry.nodes:
                if node is not None:
                    node[axis_idx] = 2 * location - node[axis_idx]
            for element in geometry.elements:
                if element is not None:
                    element.mirror_ids()
        else:
            node_map = {}
            for idx, node in enumerate(geometry.nodes):
                if node is not None:
                    # check if its close, if so, map to itself
                    if np.abs(node[axis_idx] - location) < 1e-6:
                        node_map[idx] = idx
                    else:
                        new_coords = node.copy()
                        new_coords[axis_idx] = 2 * location - new_coords[axis_idx]
                        nid = geometry.add_node(x=new_coords[0], y=new_coords[1], z=new_coords[2])
                        node_map[idx] = nid

            for elem_id, element in enumerate(geometry.elements):
                if element is not None:
                    new_node_ids = [node_map[nid] for nid in element.node_ids]
                    geometry.elements[elem_id] = element.__class__(element.element_id, new_node_ids)
                    geometry.elements[elem_id].mirror_ids()

            # copy sets
            for setname in list(geometry.node_sets.keys()):
                # create copied set
                new_set_name = setname + "_mirrored"
                geometry.add_node_set(new_set_name)
                for nid in geometry.node_sets[setname]:
                    geometry.add_node_to_set(new_set_name, node_map[nid])

            for setname in list(geometry.elem_sets.keys()):
                # create copied set
                new_set_name = setname + "_mirrored"
                geometry.add_element_set(new_set_name)
                for eid in geometry.elem_sets[setname]:
                    geometry.add_element_to_set(new_set_name, eid)

        return geometry

    def plot_2d(self):
        if self.dimension != 2:
            raise ValueError("Plotting is only supported for 2D geometries.")

        for elem in self.elements:
            if elem is not None:
                for edge in elem.connectivity():
                    node1 = self.nodes[edge[0]]
                    node2 = self.nodes[edge[1]]
                    plt.plot([node1[0], node2[0]], [node1[1], node2[1]], 'b-')

        # aspect ratio = 1
        plt.gca().set_aspect('equal', adjustable='box')

        plt.show()

    def get_surface_nodes(self) -> list[int]:
        """
        Liefert IDs aller Knoten auf der Außenfläche, ermittelt rein über
        'lonely faces' (Faces, die nur in genau einem 3D-Solid vorkommen).
        Keine Filterung, kein candidate_set, keine Fallbacks.
        Shell-Elemente tragen ihre Knoten direkt zur Oberfläche bei.
        """
        surf_nodes = set()

        # Hilfsfunktionen nur für die 3D-Solids
        def _corner_nodes(elem):
            ids = elem.node_ids
            t = elem.elem_type
            if t in ("C3D8", "C3D20"):  # Hexaeder
                return ids[:8]
            if t in ("C3D4", "C3D10"):  # Tetraeder
                return ids[:4]
            if t in ("C3D6", "C3D15"):  # Wedge/Prisma
                return ids[:6]
            return ids

        def _faces_of(elem):
            t  = elem.elem_type
            cn = _corner_nodes(elem)
            if t in ("C3D8", "C3D20"):
                return [
                    (cn[0], cn[1], cn[2], cn[3]),
                    (cn[4], cn[5], cn[6], cn[7]),
                    (cn[0], cn[1], cn[5], cn[4]),
                    (cn[1], cn[2], cn[6], cn[5]),
                    (cn[2], cn[3], cn[7], cn[6]),
                    (cn[3], cn[0], cn[4], cn[7]),
                ]
            if t in ("C3D4", "C3D10"):
                return [
                    (cn[0], cn[1], cn[2]),
                    (cn[0], cn[1], cn[3]),
                    (cn[1], cn[2], cn[3]),
                    (cn[0], cn[2], cn[3]),
                ]
            if t in ("C3D6", "C3D15"):
                return [
                    (cn[0], cn[1], cn[2]),            # tri
                    (cn[3], cn[4], cn[5]),            # tri
                    (cn[0], cn[1], cn[4], cn[3]),     # quad
                    (cn[1], cn[2], cn[5], cn[4]),     # quad
                    (cn[2], cn[0], cn[3], cn[5]),     # quad
                ]
            return []

        # 1) Lonely-Face-Logik für 3D-Solids
        face_counter = Counter()
        for e in self.elements:
            if e is None:
                continue
            if e.elem_type.startswith("C3D"):
                for f in _faces_of(e):
                    face_counter[tuple(sorted(f))] += 1

        for key, cnt in face_counter.items():
            if cnt == 1:
                surf_nodes.update(key)

        # 2) Shells als Oberfläche (direkt alle Knoten)
        for e in self.elements:
            if e is None:
                continue
            if e.elem_type in ("S3", "S4", "S6", "S8"):
                surf_nodes.update(e.node_ids)

        return sorted(surf_nodes)

    def __str__(self):
        ret_str = ""
        ret_str += f"Dimension: {'2D' if self.dimension == 2 else '3D'}\n"
        ret_str += "Nodes    : {}\n".format(len(self.nodes))
        ret_str += "Elements : {}\n".format(len(self.elements))
        ret_str += "Node Sets:\n"
        for name, ids in self.node_sets.items():
            ret_str += "  {}: {} nodes\n".format(name, len(ids))
        ret_str += "Element Sets:\n"
        for name, ids in self.elem_sets.items():
            ret_str += "  {}: {} elements\n".format(name, len(ids))
        ret_str += "Coordinate Systems:\n"
        for name, cs in self.coordinate_systems.items():
            ret_str += f"  {name}: {cs['type']} ({len(cs['values'])} values)\n"
        ret_str += "Couplings:\n"
        for c in self.couplings:
            ret_str += ("  {type}: master={master} slave={slave} "
                        "dofs=({cx},{cy},{cz},{crx},{cry},{crz})\n").format(**c)
        ret_str += "Connectors:\n"
        for c in self.connectors:
            ret_str += f"  {c['type']}: {c['nset1']} <-> {c['nset2']} (csys={c['coord_sys']})\n"
        return ret_str

    # Re-exports / Bindings
    connectivity_node_to_element    = connectivity_node_to_element
    connectivity_element_to_element = connectivity_element_to_element
    element_element_distance_matrix = element_element_distance_matrix

    subdivided       = subdivided_geometry
    extruded         = extruded_geometry
    write_input_deck = write_input_deck
    read_input_deck  = staticmethod(read_input_deck)
    mesh_interior    = staticmethod(mesh_interior)
