from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
from scipy.sparse import lil_matrix

import numpy as np

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
