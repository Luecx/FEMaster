

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

        self.materials = {}
        self.sections  = []
        self.profiles  = []

        self.loads     = {} # set, (x,y,z,rx,ry,rz)
        self.supps     = {} # set, (x,y,z,rx,ry,rz)

        self.steps     = [] # name -> loads, supps

    def add_node(self, node_id=-1, x=0, y=0, z=0):
        # Resize nodes list if necessary
        if node_id == -1:
            node_id = len(self.nodes)
        if node_id >= len(self.nodes):
            self.nodes.extend([None] * (node_id - len(self.nodes) + 1))
        self.nodes[node_id] = (x, y, z)
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

    def add_section_shell(self, material, elset, thickness):
        self.sections.append({'type': 'SHELL', 'material': material, 'elset': elset, 'thickness': thickness})

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

    def add_load(self, name, set, fx=None, fy=None, fz=None, mx=None, my=None, mz=None):
        self.loads[name] = {'set':set, 'data': (fx, fy, fz, mx, my, mz)}

    def add_supp(self, name, set, fx=None, fy=None, fz=None, mx=None, my=None, mz=None):
        self.supps[name] = {'set':set, 'data': (fx, fy, fz, mx, my, mz)}

    def add_step(self, type, name="Step-1", **kwargs):
        step = {'name': name, 'type': type.upper()}

        if step['type'] == "LINEAR STATIC":
            required_keys = ['loads', 'supps']
        elif step['type'] == "EIGENFREQ":
            required_keys = ['supps', 'numeigenvalues']
        else:
            raise ValueError(f"Unsupported step type: {type}")

        for key in required_keys:
            if key not in kwargs:
                raise ValueError(f"{key} must be provided for step type '{type}'")

        step.update(kwargs)
        self.steps.append(step)


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
        geometry.nodes = self.nodes.copy()
        geometry.node_sets = {name: ids.copy() for name, ids in self.node_sets.items()}
        geometry.elem_sets = {name: ids.copy() for name, ids in self.elem_sets.items()}
        geometry.elements = copy.deepcopy(self.elements)

        edge_midpoints = geometry.compute_edge_midpoints()

        # Convert elements to second order and replace inplace
        for i in range(len(self.elements)):
            element = self.elements[i]
            if element is not None:
                geometry.elements[i] = element.to_second_order(edge_midpoints)
        return geometry

    def mirrored(self, plane, location, copy=True):
        geometry = Geometry(self.dimension)
        geometry.nodes = self.nodes.copy()
        geometry.node_sets = {name: ids.copy() for name, ids in self.node_sets.items()}
        geometry.elem_sets = {name: ids.copy() for name, ids in self.elem_sets.items()}
        geometry.elements = copy.deepcopy(self.elements)

        axis_idx = ['yz', 'xz', 'xy'].index(plane)

        if not copy:
            # mirror all nodes and elements inplace
            for node in geometry.nodes:
                if node is not None:
                    node[axis_idx] = 2 * location - node[axis_idx]
            for element in geometry.elements:
                if element is not None:
                    element.mirror_ids()


        else:
            node_map = {}
            for idx,node in enumerate(geometry.nodes):
                if node is not None:

                    # check if its close, if so, map to itself
                    if np.abs(node[axis_idx] - location) < 1e-6:
                        node_map[idx] = idx
                    else:
                        new_coords = node.copy()
                        new_coords[axis_idx] = 2 * location - new_coords[axis_idx]
                        nid = geometry.add_node(x = new_coords[0], y = new_coords[1], z = new_coords[2])
                        node_map[idx] = nid

            for elem_id, element in enumerate(geometry.elements):
                if element is not None:
                    new_node_ids = [node_map[nid] for nid in element.node_ids]
                    geometry.elements[elem_id] = element.__class__(element.element_id, new_node_ids)
                    geometry.elements[elem_id].mirror_ids()

            # copy sets
            for setname in geometry.node_sets:
                # create copied set
                new_set_name = setname + "_mirrored"
                geometry.add_node_set(new_set_name)
                for idx, nid in enumerate(geometry.node_sets[setname]):
                    geometry.add_node_to_set(new_set_name, node_map[nid])

            for setname in geometry.elem_sets:
                # create copied set
                new_set_name = setname + "_mirrored"
                geometry.add_element_set(new_set_name)
                for idx, eid in enumerate(geometry.elem_sets[setname]):
                    geometry.add_element_to_set(new_set_name, eid)


        return geometry


        # # check if axis is correct
        # if plane not in ['yz', 'xz', 'xy']:
        #     raise ValueError(f"Unsupported plane: {plane}, use 'yz', 'xz' or 'xy'")
        #
        # # copy all nodes that are outside of the plane
        # for node_id in range(len(self.nodes)):
        #     node = self.nodes[node_id]
        #     if node is None:
        #         continue
        #
        #     if node[axis_idx] - location != 0:
        #         geometry.add_node(node_id, *node)

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
        return ret_str

    connectivity_node_to_element = connectivity_node_to_element
    connectivity_element_to_element = connectivity_element_to_element
    element_element_distance_matrix = element_element_distance_matrix


    subdivided          = subdivided_geometry
    extruded            = extruded_geometry
    write_input_deck    = write_input_deck
    read_input_deck     = staticmethod(read_input_deck)
    mesh_interior       = staticmethod(mesh_interior)

