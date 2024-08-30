import numpy as np
import matplotlib.pyplot as plt

class SpringSimulation:
    def __init__(self, nodes, edges):

        # make sure nodes has no z component
        if len(nodes[0]) == 3:
            nodes = [[x, y] for x, y, z in nodes]

        self.nodes = np.array(nodes)
        self.edges = edges
        self.boundary_nodes = self._initialize_boundary_info()
        self.edges          = self._unique_edges(edges)
        self.adjacency_list = self._build_adjacency()
        self._sort_adjacency()

        self.divisions       = self._compute_divisions()
        self.spanning_angles = self._compute_spanned_angles()
        self.target_angles   = self._compute_target_angles()


    def _build_adjacency(self):
        # Adjacency list for each node containing a list of all connected nodes
        adjacency = [[] for _ in range(len(self.nodes))]
        for edge in self.edges:
            adjacency[edge[0]].append(edge[1])
            adjacency[edge[1]].append(edge[0])
        return adjacency
    def _unique_edges(self, edges):
        unique_edges = set()
        for edge in edges:
            edge = tuple(sorted(edge))
            unique_edges.add(edge)
        return list(unique_edges)
    def _initialize_boundary_info(self):
        boundary_nodes = np.zeros(len(self.nodes), dtype=bool)

        edge_occurrences = {}

        for edge in self.edges:
            edge = tuple(sorted(edge))
            if edge in edge_occurrences:
                edge_occurrences[edge] += 1
            else:
                edge_occurrences[edge] = 1

        for edge, occurrences in edge_occurrences.items():
            if occurrences == 1:
                boundary_nodes[edge[0]] = True
                boundary_nodes[edge[1]] = True

        self.edge_occurrences = edge_occurrences

        return boundary_nodes
    def _sort_adjacency(self):
        for node, neighbors in enumerate(self.adjacency_list):
            if len(neighbors) > 1:
                center = self.nodes[node]
                coordinates = self.nodes[neighbors]
                angles = np.arctan2(coordinates[:, 1] - center[1], coordinates[:, 0] - center[0])
                sorted_neighbors = [x for _, x in sorted(zip(angles, neighbors))]

                if self.boundary_nodes[node]:
                    boundary_neighbors = [n for n in sorted_neighbors if self.boundary_nodes[n]]
                    # find two adjacent boundary nodes and start the list with the second one
                    for n in range(len(sorted_neighbors)):
                        if n in boundary_neighbors and (n + 1) % len(boundary_neighbors) in boundary_neighbors:
                            sorted_neighbors = sorted_neighbors[n:] + sorted_neighbors[:n]
                            break

                    # if len(boundary_neighbors) == 2:
                    #     if self.is_boundary_node(sorted_neighbors[0]) and self.is_boundary_node(sorted_neighbors[-1]):
                    #         pass
                    #     else:
                    #         start_index = sorted_neighbors.index(boundary_neighbors[1])
                    #         sorted_neighbors = sorted_neighbors[start_index:] + sorted_neighbors[:start_index]
                    # else:
                    #     print(node)
                    #     print(self.nodes[node])
                    #     print(boundary_neighbors)
                    #     pass
                    #     # raise ValueError("Boundary node has something other than 2 boundary neighbors")
                self.adjacency_list[node] = sorted_neighbors
    def _compute_spanned_angles(self):
        # for each node, compute the total spanned angle by its neighbors. its the sum of compute_enclosing_angles
        spanning_angles = [[] for _ in range(len(self.nodes))]
        enclosed_angles = self.compute_enclosing_angles()
        for node, neighbors in enumerate(self.adjacency_list):
            if len(neighbors) > 1:
                spanning_angles[node] = sum(enclosed_angles[node])
        return spanning_angles
    def _compute_divisions(self):
        # how many divisions, there are. for boundary nodes, it is neighbours - 1, for others it is neighbours
        divisions = [[] for _ in range(len(self.nodes))]
        for node, neighbors in enumerate(self.adjacency_list):
            if len(neighbors) > 1:
                divisions[node] = len(neighbors) - 1 if self.is_boundary_node(node) else len(neighbors)
        return divisions
    def _compute_target_angles(self):
        # target angles = spanning angles / divisions
        target_angles = [[] for _ in range(len(self.nodes))]
        spanning_angles = self.spanning_angles
        divisions = self.divisions
        for node in range(len(self.nodes)):
            target_angles[node] = spanning_angles[node] / divisions[node]
        return target_angles



    def is_boundary_node(self, node):
        return self.boundary_nodes[node]

    def compute_edge_distances(self):
        return {edge: np.linalg.norm(self.nodes[edge[0]] - self.nodes[edge[1]]) for edge in self.edges}

    def compute_mean_distance(self, node):
        neighbors = self.adjacency_list[node]
        distances = [np.linalg.norm(self.nodes[node] - self.nodes[nb]) for nb in neighbors]
        return np.mean(distances) if distances else 0

    def compute_edge_directions(self):
        return {edge: self.nodes[edge[1]] - self.nodes[edge[0]] for edge in self.edges}

    def compute_enclosing_angles(self):
        angles = [[] for _ in range(len(self.nodes))]
        for node, neighbors in enumerate(self.adjacency_list):
            if len(neighbors) > 1:
                center = self.nodes[node]
                coordinates = self.nodes[neighbors]
                neighbor_angles = np.arctan2(coordinates[:, 1] - center[1], coordinates[:, 0] - center[0])


                if not self.is_boundary_node(node):
                    neighbor_angles = np.concatenate([neighbor_angles, neighbor_angles[0:1]])
                angles[node] = np.diff(neighbor_angles)

                # any angle less than 0 should be added 2pi
                for i in range(len(angles[node])):
                    if angles[node][i] < 0:
                        angles[node][i] += 2 * np.pi

        return angles


    def create_nodal_forces_rotational_spring(self):
        angles = self.compute_enclosing_angles()
        edge_angle_adjustments = { edge : 0 for edge in self.edges }

        # for each node, go through all the edges and compute a moment
        for n in range(len(self.nodes)):

            for idx, angle in enumerate(angles[n]):

                diff_to_target =self.target_angles[n] - angle

                # strong effect to avoid too sharp angles
                if angle < 2 * np.pi * 10 / 360:
                    diff_to_target *= 3


                adj1 = self.adjacency_list[n][idx]
                adj2 = self.adjacency_list[n][(idx + 1) % len(self.adjacency_list[n])]

                edge_angle_adjustments[tuple(sorted((n, adj1)))] += diff_to_target / 2
                edge_angle_adjustments[tuple(sorted((n, adj2)))] -= diff_to_target / 2


        # compute nodal forces from angle adjustments per edge
        nodal_forces = np.zeros_like(self.nodes)
        for edge, angle_adjustment in edge_angle_adjustments.items():
            edge_center = (self.nodes[edge[0]] + self.nodes[edge[1]]) / 2
            edge_length = np.linalg.norm(self.nodes[edge[0]] - self.nodes[edge[1]])

            # node 1 adjust
            node_1_vector = self.nodes[edge[0]] - edge_center
            node_2_vector = self.nodes[edge[1]] - edge_center

            # compute the direction of movement (perpendicular to the edge)
            node_1_direction = np.array([node_1_vector[1], -node_1_vector[0]])
            node_2_direction = np.array([node_2_vector[1], -node_2_vector[0]])

            # adjust the position by direction * angle_adjustment * edge_length
            nodal_forces[edge[0]] += node_1_direction * angle_adjustment
            nodal_forces[edge[1]] += node_2_direction * angle_adjustment
        return nodal_forces

    def create_nodal_forces_directed_spring(self):
        edge_length_changes = {edge: 0 for edge in self.edges}
        for node in range(len(self.nodes)):

            if not self.is_boundary_node(node):
                neighbors = self.adjacency_list[node]
                directions = self.nodes[neighbors] - self.nodes[node]
                distances = np.linalg.norm(directions, axis=1)

                target_distance = 0

                for idx, neighbor in enumerate(neighbors):
                    edge = tuple(sorted((node, neighbor)))
                    edge_length_changes[edge] += (distances[idx] - target_distance)

        nodal_forces = np.zeros_like(self.nodes)
        # adjust the node positions
        for edge, length_change in edge_length_changes.items():
            edge_vector = self.nodes[edge[1]] - self.nodes[edge[0]]
            edge_direction = edge_vector / np.linalg.norm(edge_vector)
            nodal_forces[edge[0]] += edge_direction * length_change
            nodal_forces[edge[1]] -= edge_direction * length_change

        return nodal_forces

    def create_nodal_forces(self, rotational_spring_factor=1.0, directed_spring_factor=1.0):
        return self.create_nodal_forces_directed_spring() * directed_spring_factor \
             + self.create_nodal_forces_rotational_spring() * rotational_spring_factor

    def step(self, rotation_factor=0.1, length_factor=0.1):
        forces = self.create_nodal_forces(rotation_factor, length_factor)

        # adjust all non boundary nodes by the forces
        current_pos = self.nodes.copy()
        for node in range(len(self.nodes)):
            if not self.is_boundary_node(node):
                self.nodes[node] += forces[node]
        difference = np.linalg.norm(self.nodes - current_pos)
        return difference / len(self.nodes)

    def plot(self, show_movement=True):

        if show_movement:
            forces = self.create_nodal_forces()

        for edge in self.edges:
            # if the edge occured only once, make it blue, if it occured twice, make it black, if it occured more than twice, make it red
            color = 'blue' if self.edge_occurrences[tuple(sorted(edge))] == 1 else 'black' if self.edge_occurrences[tuple(sorted(edge))] == 2 else 'red'

            plt.plot([self.nodes[edge[0]][0], self.nodes[edge[1]][0]],
                     [self.nodes[edge[0]][1], self.nodes[edge[1]][1]], color)
        # plot all nodes
        plt.scatter(self.nodes[:, 0], self.nodes[:, 1])

        # plot movement of nodes as arrows
        if show_movement:
            for node in range(len(self.nodes)):
                if not self.is_boundary_node(node):
                    plt.arrow(self.nodes[node][0], self.nodes[node][1], forces[node][0], forces[node][1], color='red', head_width=0.05)

        # plot boundary nodes with another color
        boundary_nodes = self.nodes[self.boundary_nodes]
        plt.scatter(boundary_nodes[:, 0], boundary_nodes[:, 1], color='red')
        # aspect ratio = 1
        plt.gca().set_aspect('equal', adjustable='box')
        plt.show()

    def __str__(self):
        return f"SpringSimulation with {len(self.nodes)} nodes and {len(self.edges)} edges."

# nodes = [[0, 0], [1, 0], [1, 1], [0, 1], [0.6, 0.5]]
#
# edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 3), (0, 4), (1, 4), (0, 4), (1, 4), (4, 3), (4,2), (2,4)]
#
# simulation = SpringSimulation(nodes, edges)
#
# simulation.plot()
# for _ in range(10):
#     simulation.step()
# simulation.plot()

# Example usage:
# Assuming nodes and edges are defined:
# nodes = [[x, y], ...]
# edges = [(0, 1), (1, 2), ...]
# simulation = SpringSimulation(nodes, edges)
# simulation.adjust_node_positions(10, damping=0.1)
