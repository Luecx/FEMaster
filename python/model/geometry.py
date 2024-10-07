from .elements import *

from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix

import re
import numpy as np
import copy
import gmsh
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
    }

    def __init__(self, dimension=3):
        self.dimension = dimension  # 2 for 2D, 3 for 3D
        self.nodes = []
        self.elements = []
        self.node_sets = {"NALL": []}
        self.elem_sets = {"EALL": []}

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

    def determine_element_order(self):
        for element in self.elements:
            if element is not None:
                if element.elem_type in ['C2D3', 'C2D4', 'C3D4', 'C3D6', 'C3D8']:
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

    def get_boundary_node_ids(self):
        if self.dimension != 2:
            raise ValueError("Boundary node IDs can only be determined for 2D geometries.")

        edges = self.get_edges()
        boundary_nodes = set()

        adjacency = {node : [] for node in range(len(self.nodes))}
        for edge in edges:
            adjacency[edge[0]].append(edge[1])
            adjacency[edge[1]].append(edge[0])

        for i in range(len(self.nodes)):
            if self.nodes[i] is None:
                continue

            adjacent_nodes = adjacency[i]
            adjacency_set = {}

            for node in adjacent_nodes:
                their_adjas = adjacency[node]
                for n in their_adjas:
                    if n in adjacent_nodes:
                        if n not in adjacency_set:
                            adjacency_set[n] = 1
                        else:
                            adjacency_set[n] += 1

            for node in adjacent_nodes:
                if adjacency_set[node] != 2:
                    boundary_nodes.add(i)
        return boundary_nodes

    def equalize_node_spacing(self, iterations, strength = 0.1, min_boundary_distance = 0.1):
        edges = self.get_edges()

        adjacency = {node : [] for node in range(len(self.nodes))}
        node_surrounding_egdes = {node : [] for node in range(len(self.nodes))}

        for edge in edges:
            adjacency[edge[0]].append(edge[1])
            adjacency[edge[1]].append(edge[0])

        for node in adjacency:
            # sort the nodes clockwise
            node_positions = [np.array(self.nodes[adj_node]) for adj_node in adjacency[node]]
            node_positions = np.array(node_positions)
            center = np.mean(node_positions, axis=0)
            angles = np.arctan2(node_positions[:, 1] - center[1], node_positions[:, 0] - center[0])
            adjacency[node] = [adj_node for _, adj_node in sorted(zip(angles, adjacency[node]))]

            # fill in node surrounding edges
            for i in range(len(adjacency[node])):
                n1 = adjacency[node][i]
                n2 = adjacency[node][(i + 1) % len(adjacency[node])]
                if (n1, n2) in edges:
                    node_surrounding_egdes[node].append(tuple(sorted([n1, n2])))
                if (n2, n1) in edges:
                    node_surrounding_egdes[node].append(tuple(sorted([n2, n1])))

        boundary_nodes = self.get_boundary_node_ids()
        print("found ", len(boundary_nodes), " boundary nodes")

        # boolean mask for boundary nodes to index by node id
        boundary_nodes = [node in boundary_nodes for node in range(len(self.nodes))]

        for _ in range(iterations):
            edge_forces = {edge : 0 for edge in edges}

            edge_directions = {edge : (np.array(self.nodes[edge[1]]) - np.array(self.nodes[edge[0]])) for edge in edges}
            edge_distances  = {edge : np.linalg.norm(edge_directions[edge]) for edge in edges}

            for node in range(len(self.nodes)):

                adjacent_ids = adjacency[node]

                # theoretical_center = np.mean([np.array(self.nodes[adj_node]) for adj_node in adjacent_ids], axis=0)
                distances = [edge_distances[tuple(sorted([node, adj_node]))] for adj_node in adjacent_ids]

                # distances = [edge_distances[tuple(sorted([node, adj_node]))] for adj_node in adjacent_ids]
                mean_distance = np.mean(distances)

                # assign forces based on the difference between the mean distance and the actual distance
                for idx, adj_node in enumerate(adjacent_ids):
                    edge = tuple(sorted([node, adj_node]))
                    edge_forces[edge] += (mean_distance - distances[idx]) / 2

            new_node_pos = self.nodes.copy()

            # adjust positions, go through each edge
            for edge in edges:
                dir   = edge_directions[edge]
                dist  = edge_distances[edge]
                force = edge_forces[edge] * strength

                n1, n2 = edge
                movement = force * dir / dist


                if not boundary_nodes[n1]:

                    # new_pos = np.array(new_node_pos[n1]) - movement

                    # for edge in node_surrounding_egdes[n1]:
                    #     #

                    new_node_pos[n1] = tuple(np.array(new_node_pos[n1]) - movement)
                if not boundary_nodes[n2]:
                    new_node_pos[n2] = tuple(np.array(new_node_pos[n2]) + movement)

            movement_change = 0
            for i in range(len(self.nodes)):
                movement_change += np.linalg.norm(np.array(self.nodes[i]) - np.array(new_node_pos[i]))
                movement_change /= len(self.nodes)

            print("iteration: ", _, "movement change: ", movement_change)

            self.nodes = new_node_pos

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

    def extruded(self, n, spacing=1):
        # Determine element order
        element_order = self.determine_element_order()
        if element_order == 0:
            raise ValueError("Element order could not be determined or no elements are present.")

        geometry = Geometry(3)  # New 3D geometry
        max_node_id = len(self.nodes)

        # Copy nodes and node sets
        for i in range(n * element_order + 1):
            for node_id, node in enumerate(self.nodes):
                if node is not None:
                    x, y, z = node
                    new_node_id = i * max_node_id + node_id
                    geometry.add_node(new_node_id, x, y, z + i * spacing / (2 if element_order == 2 else 1))

        new_element_id = 1
        # Copy elements and element sets
        for i in range(n):
            for element in self.elements:
                if element is not None:
                    if element_order == 1:
                        node_ids_layer1 = [i * max_node_id + node_id for node_id in element.node_ids]
                        node_ids_layer2 = [(i + 1) * max_node_id + node_id for node_id in element.node_ids]

                        new_ids = node_ids_layer1 + node_ids_layer2

                    else:
                        node_ids_layer1 = [i * max_node_id + node_id for node_id in element.node_ids]
                        node_ids_layer2 = [(i + 2) * max_node_id + node_id for node_id in element.node_ids]
                        node_ids_layerm = [(i + 1) * max_node_id + node_id for node_id in element.node_ids]

                        # remove the second half of the mid layer
                        node_ids_layerm = node_ids_layerm[:len(node_ids_layerm) // 2]

                        # split the layer1 into two parts of equal size
                        node_ids_layer1a = node_ids_layer1[:len(node_ids_layer1) // 2]
                        node_ids_layer1b = node_ids_layer1[len(node_ids_layer1) // 2:]

                        # split the layer2 into two parts
                        node_ids_layer2a = node_ids_layer2[:len(node_ids_layer2) // 2]
                        node_ids_layer2b = node_ids_layer2[len(node_ids_layer2) // 2:]

                        new_ids = node_ids_layer1a + node_ids_layer2a + node_ids_layer1b + node_ids_layer2b + node_ids_layerm


                    geometry.add_element(new_element_id, "C3D"+str(len(new_ids)), new_ids)
                    new_element_id += 1
        # Copy sets
        for name, ids in self.node_sets.items():
            if name == "NALL":
                continue
            geometry.add_node_set(name)
            for node_id in ids:
                for i in range(n + 1):
                    geometry.node_sets[name].append(i * max_node_id + node_id)

        for name, ids in self.elem_sets.items():
            if name == "EALL":
                continue
            geometry.add_element_set(name)
            for elem_id in ids:
                for i in range(n):
                    geometry.elem_sets[name].append(i * len(self.elements) + elem_id)

        return geometry

    def subdivided(self, n=1, only_quads=False):
        # Determine element order
        element_order = self.determine_element_order()
        if element_order == 0:
            raise ValueError("Element order could not be determined or no elements are present.")
        if element_order == 2:
            raise ValueError("Subdivision is not supported for second-order elements.")

        geometry = Geometry(self.dimension)
        geometry.nodes = self.nodes.copy()
        geometry.node_sets = {name: ids.copy() for name, ids in self.node_sets.items()}
        geometry.elem_sets = {name: ids.copy() for name, ids in self.elem_sets.items()}
        geometry.elements = copy.deepcopy(self.elements)

        for i in range(n):
            next_geom = Geometry(geometry.dimension)
            next_geom.nodes = geometry.nodes.copy()
            next_geom.node_sets = {name: ids.copy() for name, ids in geometry.node_sets.items()}
            next_geom.elem_sets = {name: ids.copy() for name, ids in geometry.elem_sets.items()}
            next_geom.elements = copy.deepcopy(geometry.elements)

            midpoints = next_geom.compute_edge_midpoints()

            for i in range(len(next_geom.elements)):
                element = next_geom.elements[i]
                if element is not None:
                    element.subdivide(midpoints, next_geom, only_quads=only_quads)

            geometry = next_geom
        return geometry

    def write_input_deck(self, filename):
        with open(filename, 'w') as file:
            # Write nodes
            file.write("*NODE, NSET=NALL\n")
            for id, node in enumerate(self.nodes):
                if node is not None:
                    file.write(f"{id},{node[0]},{node[1]},{node[2]}\n")

            # Write elements by type
            element_classes = {k : [] for k in Geometry.element_classes.keys()}

            for element in self.elements:
                if element is not None:
                    element_classes[element.elem_type].append(element)

            for elem_type, elements in element_classes.items():
                if elements:
                    file.write(f"*ELEMENT, TYPE={elem_type}, ELSET=EALL\n")
                    for element in elements:
                        nodes = ','.join(map(str, [n for n in element.node_ids]))
                        file.write(f"{element.element_id},{nodes}\n")

            # Write node sets
            for name, ids in self.node_sets.items():
                if name != "NALL":
                    file.write(f"*NSET, NAME={name}\n")
                    file.write(',\n'.join(map(str, [id for id in ids])) + '\n')

            # Write element sets
            for name, ids in self.elem_sets.items():
                if name != "EALL":
                    file.write(f"*ELSET, NAME={name}\n")
                    file.write(',\n'.join(map(str, [id for id in ids])) + '\n')

    @staticmethod
    def read_input_deck(filename):
        geometry = Geometry()
        with open(filename, 'r') as file:
            lines = file.readlines()

        def preprocess_line(line):
            """Convert the line to upper case, remove spaces, and replace them with commas."""
            return re.sub(r'\s+', ',', line.strip().upper())

        def is_command(line):
            """Check if the line starts with '*' indicating a command."""
            return line.startswith('*')

        def require(key_names, parts, default=None, options=None):
            """Extract a key value from the parts based on the provided key names."""
            for key in key_names:
                if key in parts:
                    key_index = parts.index(key)
                    if key_index + 1 < len(parts):
                        if options is not None and parts[key_index + 1] not in options:
                            raise ValueError(f"Invalid value for {key}: {parts[key_index + 1]} \n\tExpected: {options}")
                        return parts[key_index + 1]
            return default

        def parse_line(line, dtype=float):
            """Parse a line of data into a list of values, converting to a specified type."""
            line = re.sub(r',\s*$', '', line)  # Remove trailing commas and spaces
            return list(map(dtype, re.split(r'[,\s]+', line.strip())))

        key_word        = None
        current_nset    = None
        current_elset   = None
        index = 0

        while index < len(lines):
            line = preprocess_line(lines[index])

            if not line or line.startswith('**'):  # Skip comments or empty lines
                index += 1
                continue

            if is_command(line):
                parts = re.split(r'[,=\s]+', line)
                key_word = parts[0].upper()

                # Process different commands
                if key_word == '*NODE':
                    current_nset = require(["NSET", "NAME"], parts, default='NALL')
                    geometry.add_node_set(current_nset)

                elif key_word == '*ELEMENT':
                    elem_type = require(["TYPE"], parts, options=Geometry.element_classes.keys())
                    current_elset = require(["ELSET", "NAME"], parts, default='EALL')
                    geometry.add_element_set(current_elset)

                elif key_word == '*NSET':
                    current_nset = require(["NAME", "NSET"], parts)
                    geometry.add_node_set(current_nset)

                elif key_word == '*ELSET':
                    current_elset = require(["NAME", "ELSET"], parts)
                    geometry.add_element_set(current_elset)

                else:
                    key_word = None  # Unknown command
                index += 1

            else:
                if key_word is None:
                    index += 1
                    continue

                data = parse_line(line)

                # Process node data
                if key_word == '*NODE':
                    node_id, *coords = data
                    geometry.add_node(int(node_id), *coords)
                    if current_nset:
                        geometry.add_node_to_set(current_nset, int(node_id))

                # Process element data
                elif key_word == '*ELEMENT' and elem_type in Geometry.element_classes:
                    required_nodes = Geometry.element_classes[elem_type].num_nodes
                    element_data = data

                    # Read more lines if the number of nodes is not enough
                    while len(element_data) - 1 < required_nodes:
                        index += 1
                        extra_data = parse_line(preprocess_line(lines[index]))
                        element_data.extend(extra_data)

                    element_id, *node_ids = map(int, element_data)
                    geometry.add_element(int(element_id), elem_type, [nid for nid in node_ids])
                    if current_elset:
                        geometry.add_element_to_set(current_elset, int(element_id))

                # Process NSET data
                elif key_word == '*NSET' and current_nset:
                    for nid in data:
                        geometry.add_node_to_set(current_nset, int(nid))

                # Process ELSET data
                elif key_word == '*ELSET' and current_elset:
                    for eid in data:
                        geometry.add_element_to_set(current_elset, int(eid))

                index += 1

        return geometry

    def connectivity_node_to_element(self):
        node_elements = {i: [] for i in range(len(self.nodes))}

        for i, element in enumerate(self.elements):
            if element is not None:
                for node_id in element.node_ids:
                    node_elements[node_id].append(i)
        return node_elements

    def connectivity_element_to_element(self):
        node_elements = self.connectivity_node_to_element()
        n_elements = len(self.elements)

        # Initialize a sparse adjacency matrix (lil_matrix allows efficient incremental construction)
        element_elements = lil_matrix((n_elements, n_elements), dtype=np.int8)

        for node, elements in node_elements.items():
            for element in elements:
                for connected_element in elements:
                    if element != connected_element:
                        element_elements[element, connected_element] = 1

        return element_elements.tocsr()  # Convert to CSR format for more efficient access later

    def element_element_distance_matrix(self):
        element_elements  = self.connectivity_element_to_element()
        element_midpoints = self.compute_element_midpoints()
        n_elements = len(self.elements)

        # List to store the row, column, and distance entries
        rows = []
        cols = []
        distances = []

        for i in range(n_elements):
            rows.append(i)
            cols.append(i)
            distances.append(0.0)  # Diagonal elements (self-distance)

            # Only compute distances for connected elements
            connected_elements = element_elements[i]
            for j in connected_elements:
                if j > i:  # Avoid redundant calculations by ensuring j > i
                    distance = np.linalg.norm(element_midpoints[i] - element_midpoints[j])
                    rows.append(i)
                    cols.append(j)
                    distances.append(distance)

                    # Symmetric matrix, so mirror the distance
                    rows.append(j)
                    cols.append(i)
                    distances.append(distance)

        # Create a sparse matrix in COO format and convert to CSR for efficiency
        element_distances = coo_matrix((distances, (rows, cols)), shape=(n_elements, n_elements)).tocsr()
        return element_distances

    def plot_2d(self):
        if self.dimension != 2:
            raise ValueError("Plotting is only supported for 2D geometries.")

        for elem in self.elements:
            if elem is not None:
                for edge in elem.connectivity():
                    node1 = self.nodes[edge[0]]
                    node2 = self.nodes[edge[1]]
                    plt.plot([node1[0], node2[0]], [node1[1], node2[1]], 'b-')

        # def elements_to_triangles(elements):
        #     return np.array([[e.node_ids[0], e.node_ids[1], e.node_ids[2]] for e in elements])
        #
        # nodes = np.array(self.nodes)
        # elements = np.array(elements_to_triangles(self.elements))
        #
        # # plot the geometry
        # fig, ax = plt.subplots()
        # ax.set_aspect('equal')
        # ax.set_xlim(np.min(nodes[:, 0]) - 1, np.max(nodes[:, 0]) + 1)
        # ax.set_ylim(np.min(nodes[:, 1]) - 1, np.max(nodes[:, 1]) + 1)
        #
        #
        # triangulation = tri.Triangulation(nodes[:, 0], nodes[:, 1], elements)
        # ax.triplot(triangulation, color='blue')
        plt.show()

    @staticmethod
    def mesh_interior(boundary_points, second_order=False, mesh_type=0):
        """
        Generate a 2D mesh for a given set of boundary points.

        Parameters
        ----------
        boundary_points : list of list
            List of 2D points defining the boundary of the mesh in [x, y] format.
        second_order : bool, optional
            Whether to generate second-order elements. Default is False.
        mesh_type : int, optional
            Type of elements to generate:
                0 - Triangles only (default)
                1 - Mixed mesh (triangles and quadrilaterals)
                2 - Quadrilaterals only

        Returns
        -------
        geometry : Geometry object
            Geometry object containing the generated mesh nodes and elements.

        Notes
        -----
        Gmsh will generate 9-node quadrilateral elements for second-order quads by default.
        This function automatically converts 9-node elements to 8-node elements by removing
        the center node. The resulting mesh will always contain C2D3, C2D4, C2D6, or C2D8 elements.
        """

        import gmsh
        import numpy as np

        # Remove last point if it's a duplicate of the first point
        if np.linalg.norm(boundary_points[0] - boundary_points[-1]) < 1e-6:
            boundary_points = boundary_points[:-1]

        gmsh.initialize()
        gmsh.model.add("Mesh")

        # Add points
        point_tags = []
        for pt in boundary_points:
            tag = gmsh.model.geo.addPoint(pt[0], pt[1], 0)
            point_tags.append(tag)

        # Add lines between consecutive points
        line_tags = []
        num_points = len(point_tags)
        for i in range(num_points):
            tag = gmsh.model.geo.addLine(point_tags[i], point_tags[(i + 1) % num_points])
            line_tags.append(tag)

        # Create a closed loop and surface
        curve_loop_tag = gmsh.model.geo.addCurveLoop(line_tags)
        surface_tag = gmsh.model.geo.addPlaneSurface([curve_loop_tag])

        # Synchronize before meshing
        gmsh.model.geo.synchronize()

        # Set global meshing options
        gmsh.option.setNumber("Mesh.Algorithm", 6)  # Delaunay meshing for 2D
        gmsh.option.setNumber("Mesh.ElementOrder", 1 if not second_order else 2)  # First or second order elements

        # Handle mesh type: 0 = Triangles only, 1 = Mixed (tri + quad), 2 = Quad only
        if mesh_type in {1, 2}:  # Enable recombination for mixed or quad-only mesh
            gmsh.model.geo.mesh.setRecombine(2, surface_tag)

        # Synchronize again after recombination settings
        gmsh.model.geo.synchronize()

        # Generate the mesh
        gmsh.model.mesh.generate(2)

        # Extract nodes and elements
        node_data = gmsh.model.mesh.getNodes()
        nodes = node_data[1].reshape((-1, 3))

        element_types, element_tags, node_tags_flattened = gmsh.model.mesh.getElements(dim=2)

        # Manually define the number of nodes per element type
        element_node_count = {
            2: 3,  # 3-node triangle
            3: 4,  # 4-node quadrilateral
            9: 6,  # 6-node second-order triangle
            10: 9,  # 9-node second-order quadrilateral
        }

        # Convert flattened node tags to list of lists of node tags per element
        node_tags_per_element = []
        for elem_type, node_tags in zip(element_types, node_tags_flattened):
            num_nodes_per_elem = element_node_count.get(elem_type, 0)
            if num_nodes_per_elem > 0:
                node_tags_per_element.append(np.array(node_tags).reshape(-1, num_nodes_per_elem))

        # Initialize Geometry object
        geometry = Geometry(dimension=2)

        # Map node ids to indices
        node_ids = {node_id: idx for idx, node_id in enumerate(node_data[0])}

        # Add nodes to the geometry
        for idx, node in enumerate(nodes):
            geometry.add_node(idx + 1, node[0], node[1])

        elem_id = 1
        # Add elements to the geometry
        for elem_tags, nodes_per_elem in zip(element_tags, node_tags_per_element):
            for elem_tag, node_tags in zip(elem_tags, nodes_per_elem):
                node_indices = [node_ids[n] + 1 for n in node_tags]

                # Identify element type based on the number of nodes and create geometry
                if len(node_tags) == 3:
                    element_type = 'C2D3'  # 3-node triangle
                elif len(node_tags) == 4:
                    element_type = 'C2D4'  # 4-node quadrilateral
                elif len(node_tags) == 6:
                    element_type = 'C2D6'  # 6-node second-order triangle
                elif len(node_tags) == 9:
                    element_type = 'C2D8'  # Always convert 9-node to 8-node quadrilateral
                    node_indices = node_indices[:8]
                else:
                    raise ValueError(f"Unsupported element type with {len(node_tags)} nodes.")

                geometry.add_element(elem_id, element_type, node_indices)
                elem_id += 1

        # Finalize Gmsh
        gmsh.finalize()
        return geometry


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

