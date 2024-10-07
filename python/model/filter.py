
import numpy as np
from scipy import sparse
from enum import Enum
from scipy.spatial import KDTree
from .symmetry import apply_symmetry

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
