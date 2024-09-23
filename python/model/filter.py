from .symmetry import apply_symmetry

import numpy as np
from scipy import sparse
from scipy.sparse import coo_matrix, find
from enum import Enum
from tqdm import tqdm
import time

from scipy.spatial import KDTree
# from sklearn.neighbors import KDTree


class FilterFunction(Enum):
    CONSTANT = 1
    GAUSSIAN = 2
    LINEAR   = 3

class Filter:
    def __init__(self, coords, sigma, symmetries={}, filter_func=FilterFunction.GAUSSIAN, epsilon=1e-3):
        self.sigma = sigma

        self.filter_func = filter_func
        self.symmetries = symmetries
        self.epsilon = epsilon

        # Precompute the coordinates with symmetries applied
        self.n_coords = coords.shape[0]
        self.coords = self._apply_symmetries(coords=coords)
        self.n_coords_internal = self.coords.shape[0]
        self.ghost_factor = self.n_coords_internal // self.n_coords

        # Build the KD-Tree
        self.kd_tree = KDTree(self.coords)

        # Compute the influence radius based on the filter function
        self.influence_radius = self._compute_influence_radius()

        # Precompute the sparse weight matrix
        self.weights = self._compute_weights()

    def _compute_influence_radius(self):
        if self.filter_func == FilterFunction.CONSTANT:
            radius = self.sigma
        elif self.filter_func == FilterFunction.LINEAR:
            radius = self.sigma
        elif self.filter_func == FilterFunction.GAUSSIAN:
            radius = 3 * self.sigma
        return radius

    def _apply_symmetries(self, coords):
        # Assume `apply_symmetry` is a function that applies symmetries to the coordinates
        new_coords, _ = apply_symmetry(coords=coords,
                                       values=np.zeros(self.n_coords),
                                       symmetries=self.symmetries)
        return new_coords


    def _compute_weights(self):

        if self.influence_radius == 0:
            return sparse.eye(self.n_coords_internal)

        distances = self.kd_tree.sparse_distance_matrix(self.kd_tree, self.influence_radius)
        weights   = distances.tocsr()

        bytes = weights.data.nbytes + weights.indptr.nbytes + weights.indices.nbytes
        mbytes = bytes / 1024 / 1024
        print("memory footprint of weights: ", mbytes, " MBytes")

        # apply filter to the weights
        if self.filter_func == FilterFunction.CONSTANT:
            # set all non zero values to 1
            weights.data = np.ones_like(weights.data)
        elif self.filter_func == FilterFunction.LINEAR:
            # set all non zero values to (sigma - distance) / sigma
            weights.data = np.maximum(0.0, (self.sigma - weights.data) / self.sigma)
        elif self.filter_func == FilterFunction.GAUSSIAN:
            # set all non zero values to exp(-distance^2 / (2 * sigma^2))
            weights.data = np.exp(-weights.data**2 / (2 * self.sigma**2))

        # Normalize rows to sum to 1
        row_sums = np.array(weights.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1.0  # Prevent division by zero
        scaling_factors = sparse.diags(1.0 / row_sums)
        weights = scaling_factors.dot(weights)

        min_nb = np.min(np.diff(weights.indptr))
        max_nb = np.max(np.diff(weights.indptr))
        avg_nb = np.mean(np.diff(weights.indptr))

        print(f"Min nb: {min_nb}, Max nb: {max_nb}, Avg nb: {avg_nb:.2f}")

        return weights
        # # apply filter to distance matrix
        #
        #
        # rows = []
        # cols = []
        # data = []
        #
        # # Create a tqdm progress bar
        # tqdm_bar = tqdm(total=self.n_coords, desc="Computing weights", unit="point")
        #
        # min_nb = self.n_coords
        # max_nb = 0
        # avg_nb = 0
        #
        # # For each point, find the neighbors within the influence radius
        # for i, point in enumerate(self.coords[:self.n_coords]):
        #     indices = self.kd_tree.query_ball_point(point, self.influence_radius)
        #     # indices = self.kd_tree.query_radius(self.coords[i:i+1], self.influence_radius)[0]
        #
        #     min_nb = min(min_nb, len(indices))
        #     max_nb = max(max_nb, len(indices))
        #     avg_nb = avg_nb + (len(indices) - avg_nb) / (i + 1)
        #
        #     for j in indices:
        #         dist_squared = np.linalg.norm(point - self.coords[j])**2
        #
        #         if self.filter_func == FilterFunction.CONSTANT:
        #             weight = 1.0
        #         elif self.filter_func == FilterFunction.LINEAR:
        #             weight = max(0.0, (self.sigma - np.sqrt(dist_squared)) / self.sigma)
        #         elif self.filter_func == FilterFunction.GAUSSIAN:
        #             weight = np.exp(-dist_squared / (2 * self.sigma**2))
        #
        #         rows.append(i)
        #         cols.append(j)
        #         data.append(weight)
        #
        #     # Update the progress bar description with min, max, and avg number of neighbors
        #     tqdm_bar.set_description(f"Min nb: {min_nb}, Max nb: {max_nb}, Avg nb: {avg_nb:.2f}")
        #     tqdm_bar.update(1)
        #
        # # Close the progress bar
        # tqdm_bar.close()
        #
        # # Create a sparse matrix in CSR format
        # sparse_matrix = sparse.coo_matrix((data, (rows, cols)), shape=(self.n_coords, self.n_coords_internal)).tocsr()
        #
        # # Normalize rows to sum to 1
        # row_sums = np.array(sparse_matrix.sum(axis=1)).flatten()
        # row_sums[row_sums == 0] = 1.0  # Prevent division by zero
        # scaling_factors = sparse.diags(1.0 / row_sums)
        # sparse_matrix = scaling_factors.dot(sparse_matrix)
        #
        # return sparse_matrix

    def apply(self, values):
        """
        Apply the filter to the values using the precomputed sparse weight matrix.
        """

        new_values = np.tile(values, self.ghost_factor)
        return self.weights.dot(new_values)[:self.n_coords]

    def minimal_distance(self, threshold=1e-6):
        distances, indices = self.kd_tree.query(self.coords, k=2)
        distances = distances[:, 1]
        distances[distances < threshold] = np.inf
        return np.min(distances)



        # """
        # Compute the minimal distance to the nearest neighbor for each point using the weight matrix.
        # Only considers pairs where row != col.
        # """
        # # Extract the non-zero entries of the weight matrix
        # rows, cols, data = find(self.weights)
        #
        # # Filter out entries where row == col
        # valid_mask = rows != cols
        # filtered_rows = rows[valid_mask]
        # filtered_cols = cols[valid_mask]
        # filtered_data = data[valid_mask]
        #
        # # Find the index of the maximum weight after filtering
        # idx_max = np.argmax(filtered_data)
        # max_row = filtered_rows[idx_max]
        # max_col = filtered_cols[idx_max]
        #
        # # Compute the Euclidean distance between the corresponding coordinates
        # dist = np.linalg.norm(self.coords[max_row] - self.coords[max_col])
        # return dist
