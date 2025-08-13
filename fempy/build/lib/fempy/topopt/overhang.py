import numpy as np
import math
from .plane import Plane

class Overhang:
    def __init__(self, geometry, plane: Plane, dx=None, F=4, alpha=45, element_set="EALL", max_nodes=None):
        self.geometry = geometry
        self.plane = plane
        self.alpha = alpha
        self.F = F
        self.build_direction = plane.normal / np.linalg.norm(plane.normal)

        self.element_ids = np.array(geometry.elem_sets[element_set])
        self.element_mask = np.zeros(len(geometry.elements), dtype=bool)
        self.element_mask[self.element_ids] = True
        self.elements = [geometry.elements[i] for i in self.element_ids]

        self._project_geometry_to_rst()

        self.r_min, self.r_max = np.min(self.rst_elem[:, 0]), np.max(self.rst_elem[:, 0])
        self.s_min, self.s_max = np.min(self.rst_elem[:, 1]), np.max(self.rst_elem[:, 1])
        self.t_min, self.t_max = np.min(self.rst_elem[:, 2]), np.max(self.rst_elem[:, 2])

        if max_nodes is not None:
            dx = self._find_max_dx_for_nodes(max_nodes)
        elif dx is None:
            raise ValueError("Either dx or max_nodes must be provided.")
        self.dx = dx
        self.subspacing = dx / F
        self.dt = np.tan(np.radians(alpha)) * dx

        self.r_count, self.s_count, self.rs_grid = self._rs_grid()
        self.n_layers = int(np.ceil((self.t_max - self.t_min) / self.dt))
        self.feed_matrix = self._feed_matrix()

        print(self)
        uncovered = []
        for global_idx in self.element_ids:
            r0, r1, s0, s1, t0, t1 = self._element_bounds(global_idx)
            if r1 <= r0 or s1 <= s0 or t1 <= t0:
                uncovered.append(global_idx)
        if uncovered:
            raise RuntimeError(
                f"{len(uncovered)} elements are not covered by any grid cells. "
                f"Increase resolution (smaller dx or larger max_nodes). "
                f"Example missing ID: {uncovered[0]}"
            )


    def __str__(self):
        total_nodes = self.r_count * self.s_count * self.n_layers
        return (
            f"<Overhang>\n"
            f"  dx = {self.dx:.4f}\n"
            f"  subspacing = {self.subspacing:.4f} (F = {self.F})\n"
            f"  alpha = {self.alpha}° → dt = {self.dt:.4f}\n"
            f"  r range = [{self.r_min:.3f}, {self.r_max:.3f}], count = {self.r_count}\n"
            f"  s range = [{self.s_min:.3f}, {self.s_max:.3f}], count = {self.s_count}\n"
            f"  t range = [{self.t_min:.3f}, {self.t_max:.3f}], layers = {self.n_layers}\n"
            f"  total grid nodes = {total_nodes}\n"
        )

    def _find_max_dx_for_nodes(self, max_nodes, tol=1e-4):
        def node_count_for_dx(dx):
            subspacing = dx / self.F
            dt = np.tan(np.radians(self.alpha)) * dx
            r_count = max(1, int(np.ceil((self.r_max - self.r_min) / subspacing)))
            s_count = max(1, int(np.ceil((self.s_max - self.s_min) / subspacing)))
            n_layers = int(np.ceil((self.t_max - self.t_min) / dt))
            return r_count * s_count * n_layers

        dx_low = 1e-5
        dx_high = 1e10
        best_dx = dx_low

        while dx_high - dx_low > tol:
            mid = (dx_low + dx_high) / 2
            count = node_count_for_dx(mid)
            if count <= max_nodes:
                best_dx = mid
                dx_high = mid
            else:
                dx_low = mid
        return best_dx

    def _project_geometry_to_rst(self):
        r_axis, s_axis = self.plane.orthonormal_basis()
        origin = self.plane.point
        self.rst_node = np.zeros((len(self.geometry.nodes), 3))
        for i, pos in enumerate(self.geometry.nodes):
            if pos is not None:
                rel = np.array(pos) - origin
                self.rst_node[i, 0] = np.dot(rel, r_axis)
                self.rst_node[i, 1] = np.dot(rel, s_axis)
                self.rst_node[i, 2] = np.dot(rel, self.build_direction)

        self.rst_elem = np.zeros((len(self.elements), 3))
        for i, elem in enumerate(self.elements):
            if elem is None:
                continue
            node_ids = np.unique([
                nid for group in elem.connectivity()
                for nid in (group if isinstance(group, (list, tuple)) else [group])
            ])
            coords = self.rst_node[node_ids]
            self.rst_elem[i] = np.mean(coords, axis=0)

    def _rs_grid(self):
        r_count = max(1, int(np.ceil((self.r_max - self.r_min) / self.subspacing)))
        s_count = max(1, int(np.ceil((self.s_max - self.s_min) / self.subspacing)))
        r_vals = np.linspace(self.r_min, self.r_min + r_count * self.subspacing, r_count)
        s_vals = np.linspace(self.s_min, self.s_min + s_count * self.subspacing, s_count)
        return r_count, s_count, np.meshgrid(r_vals, s_vals, indexing='ij')

    def _feed_matrix(self):
        size = 2 * self.F + 1
        radius = np.sqrt((self.F + 0.5)**2 + 0.5**2)
        weights = np.zeros((size, size))
        samples = 10
        offsets = np.linspace(-0.5 + 0.5 / samples, 0.5 - 0.5 / samples, samples)
        for j in range(size):
            for i in range(size):
                cx, cy = i - self.F, j - self.F
                count_inside = sum(
                    (cx + dx)**2 + (cy + dy)**2 <= radius**2
                    for dx in offsets for dy in offsets
                )
                weights[j, i] = count_inside / (samples ** 2)
        return weights

    def _element_bounds(self, global_idx):
        elem = self.geometry.elements[global_idx]
        if elem is None:
            return (0, 0, 0, 0, 0, 0)

        node_ids = np.unique([
            nid for group in elem.connectivity()
            for nid in (group if isinstance(group, (list, tuple)) else [group])
        ])
        coords = self.rst_node[node_ids]

        eps_r = 0.01 * self.subspacing
        eps_s = 0.01 * self.subspacing
        eps_t = 0.01 * self.dt

        r_min, r_max = np.min(coords[:, 0]) - eps_r, np.max(coords[:, 0]) + eps_r
        s_min, s_max = np.min(coords[:, 1]) - eps_s, np.max(coords[:, 1]) + eps_s
        t_min, t_max = np.min(coords[:, 2]) - eps_t, np.max(coords[:, 2]) + eps_t

        r0 = max(math.ceil((r_min - self.r_min) / self.subspacing), 0)
        r1 = min(math.floor((r_max - self.r_min) / self.subspacing) + 1, self.r_count)
        s0 = max(math.ceil((s_min - self.s_min) / self.subspacing), 0)
        s1 = min(math.floor((s_max - self.s_min) / self.subspacing) + 1, self.s_count)
        t0 = max(math.ceil((t_min - self.t_min) / self.dt), 0)
        t1 = min(math.floor((t_max - self.t_min) / self.dt) + 1, self.n_layers)

        return r0, r1, s0, s1, t0, t1

    def mesh_to_grid(self, global_values):
        grid = np.zeros((self.r_count, self.s_count, self.n_layers))
        count = np.zeros_like(grid)

        eps_r = 0.01 * self.subspacing
        eps_s = 0.01 * self.subspacing
        eps_t = 0.01 * self.dt

        for local_idx, global_idx in enumerate(self.element_ids):
            elem = self.geometry.elements[global_idx]
            if elem is None:
                continue

            node_ids = np.unique([
                nid for group in elem.connectivity()
                for nid in (group if isinstance(group, (list, tuple)) else [group])
            ])
            coords_rst = self.rst_node[node_ids]

            r_min, r_max = np.min(coords_rst[:, 0]) - eps_r, np.max(coords_rst[:, 0]) + eps_r
            s_min, s_max = np.min(coords_rst[:, 1]) - eps_s, np.max(coords_rst[:, 1]) + eps_s
            t_min, t_max = np.min(coords_rst[:, 2]) - eps_t, np.max(coords_rst[:, 2]) + eps_t

            r0 = max(math.ceil((r_min - self.r_min) / self.subspacing), 0)
            r1 = min(math.floor((r_max - self.r_min) / self.subspacing) + 1, self.r_count)
            s0 = max(math.ceil((s_min - self.s_min) / self.subspacing), 0)
            s1 = min(math.floor((s_max - self.s_min) / self.subspacing) + 1, self.s_count)
            t0 = max(math.ceil((t_min - self.t_min) / self.dt), 0)
            t1 = min(math.floor((t_max - self.t_min) / self.dt) + 1, self.n_layers)

            for r in range(r0, r1):
                for s in range(s0, s1):
                    for t in range(t0, t1):
                        grid[r, s, t] += global_values[global_idx]
                        count[r, s, t] += 1

        with np.errstate(divide='ignore', invalid='ignore'):
            return np.divide(grid, count, out=np.zeros_like(grid), where=count > 0)

    def grid_to_mesh(self, grid):
        values = np.zeros(len(self.geometry.elements))
        counts = np.zeros(len(self.geometry.elements))

        eps_r = 0.01 * self.subspacing
        eps_s = 0.01 * self.subspacing
        eps_t = 0.01 * self.dt

        for local_idx, global_idx in enumerate(self.element_ids):
            elem = self.geometry.elements[global_idx]
            if elem is None:
                continue

            node_ids = np.unique([
                nid for group in elem.connectivity()
                for nid in (group if isinstance(group, (list, tuple)) else [group])
            ])
            coords_rst = self.rst_node[node_ids]

            r_min, r_max = np.min(coords_rst[:, 0]) - eps_r, np.max(coords_rst[:, 0]) + eps_r
            s_min, s_max = np.min(coords_rst[:, 1]) - eps_s, np.max(coords_rst[:, 1]) + eps_s
            t_min, t_max = np.min(coords_rst[:, 2]) - eps_t, np.max(coords_rst[:, 2]) + eps_t

            r0 = max(math.ceil((r_min - self.r_min) / self.subspacing), 0)
            r1 = min(math.floor((r_max - self.r_min) / self.subspacing) + 1, self.r_count)
            s0 = max(math.ceil((s_min - self.s_min) / self.subspacing), 0)
            s1 = min(math.floor((s_max - self.s_min) / self.subspacing) + 1, self.s_count)
            t0 = max(math.ceil((t_min - self.t_min) / self.dt), 0)
            t1 = min(math.floor((t_max - self.t_min) / self.dt) + 1, self.n_layers)

            for r in range(r0, r1):
                for s in range(s0, s1):
                    for t in range(t0, t1):
                        values[global_idx] += grid[r, s, t]
                        counts[global_idx] += 1

        mask = counts > 0
        values[mask] /= counts[mask]
        return values

    def get_grid_coordinates_rst(self):
        r_vals = self.rs_grid[0][:, 0]
        s_vals = self.rs_grid[1][0, :]
        t_vals = [self.t_min + i * self.dt for i in range(self.n_layers)]
        r_grid, s_grid, t_grid = np.meshgrid(r_vals, s_vals, t_vals, indexing='ij')
        return r_grid, s_grid, t_grid

    def get_grid_coordinates_xyz(self):
        r_grid, s_grid, t_grid = self.get_grid_coordinates_rst()
        r_axis, s_axis = self.plane.orthonormal_basis()
        origin = self.plane.point
        xyz = (origin +
               r_grid[..., None] * r_axis +
               s_grid[..., None] * s_axis +
               t_grid[..., None] * self.build_direction)
        return xyz

    def propagate(self, densities):
        density_grid = self.mesh_to_grid(densities)
        support_grid = np.zeros_like(density_grid)
        support_grid[:, :, 0] = density_grid[:, :, 0]

        kernel = self.feed_matrix
        k = kernel.shape[0] // 2

        for l in range(1, self.n_layers):
            for i in range(self.r_count):
                for j in range(self.s_count):
                    i0 = max(i - k, 0)
                    i1 = min(i + k + 1, self.r_count)
                    j0 = max(j - k, 0)
                    j1 = min(j + k + 1, self.s_count)

                    ki0 = k - (i - i0)
                    ki1 = ki0 + (i1 - i0)
                    kj0 = k - (j - j0)
                    kj1 = kj0 + (j1 - j0)

                    supp_patch = support_grid[i0:i1, j0:j1, l - 1]
                    dens_patch = density_grid[i0:i1, j0:j1, l - 1]
                    kern_patch = kernel[ki0:ki1, kj0:kj1]

                    combined = supp_patch * dens_patch * kern_patch
                    support_grid[i, j, l] = np.max(combined)

        support = self.grid_to_mesh(support_grid)
        return support, density_grid, support_grid

    def loss(self, densities, p=1.0, lambda_=1.0):
        """
        Compute overhang loss: penalizes material in unsupported regions.

        Parameters:
            densities (np.ndarray): Element-wise density values.
            p (float): Exponent for penalization (e.g., 1 for linear, 2 for quadratic).
            lambda_ (float): Weight for overhang loss.

        Returns:
            float: Scalar loss value.
        """
        support, _, _ = self.propagate(densities)
        return - lambda_ * np.sum(support)

    def loss_derivative(self, densities, p=1.0, lambda_=1.0):
        """
        Compute gradient of overhang loss w.r.t. densities.

        Parameters:
            densities (np.ndarray): Element-wise density values.
            p (float): Same exponent as in loss.
            lambda_ (float): Weight for overhang loss.

        Returns:
            np.ndarray: Gradient (same shape as densities).
        """
        support, _, _ = self.propagate(densities)
        grad = - lambda_ * support
        return grad
