import numpy as np
import geom
import plot
import write
import subprocess
import os
import run
import pickle
from datetime import datetime
from scipy.ndimage import gaussian_filter

class Model:
    def __init__(self, width, height, length):
        self.width  = width
        self.height = height
        self.length = length
        # geometry
        self.nodes, self.elements, self.node_lookup = self._generate_model()
        self.midpoints = self._compute_element_midpoints()

        # material
        self.youngs = 210000
        self.nu = 0.3

        # filtering
        self.proximity_radius = width / 5

        # constraints
        self.supports = {}
        self.loads    = {}

        # simp exonent
        self.exponent = 3

        # solver mode
        self.cpu = True
        self.direct = True

        # setting the path
        self.path = "topo_" + datetime.now().strftime('%Y%m%d_%H-%M-%S')

        # symmetry
        self.symmetries = {}

    def get_special_node_ids(self):
        return geom.get_special_node_ids(width=self.width, height=self.height, length=self.length, node_lookup=self.node_lookup)

    def set_proximity_radius(self, radius):
        self.proximity_radius = radius

    def set_material(self, youngs, nu):
        self.youngs = youngs
        self.nu = nu

    def add_symmetry(self, plane):
        """Add symmetry along a given plane."""
        if plane in ["xy", "yz", "xz", "x", "y", "z"]:
            self.symmetries[plane] = True
        else:
            raise ValueError(f"Invalid symmetry plane {plane}. Supported planes are 'xy', 'yz', and 'xz'.")


    def enforce_symmetry(self, values):
        """Ensure the symmetry on the input values based on specified constraints."""
        data_3d = np.reshape(values, (self.width, self.height, self.length))

        for plane in self.symmetries:
            if plane == "xy":
                flipped_data = np.flip(data_3d, axis=2)
                data_3d = 0.5 * (data_3d + flipped_data)
            elif plane == "yz":
                flipped_data = np.flip(data_3d, axis=0)
                data_3d = 0.5 * (data_3d + flipped_data)
            elif plane == "xz":
                flipped_data = np.flip(data_3d, axis=1)
                data_3d = 0.5 * (data_3d + flipped_data)
            elif plane == "x":
                rot_1 = np.rot90(data_3d, axes=(1,2))
                rot_2 = np.rot90(rot_1  , axes=(1,2))
                rot_3 = np.rot90(rot_2  , axes=(1,2))
                rot_4 = np.rot90(rot_3  , axes=(1,2))
                data_3d = 0.25 * (rot_1 + rot_2 + rot_3 + rot_4)
            elif plane == "y":
                rot_1 = np.rot90(data_3d, axes=(0,2))
                rot_2 = np.rot90(rot_1  , axes=(0,2))
                rot_3 = np.rot90(rot_2  , axes=(0,2))
                rot_4 = np.rot90(rot_3  , axes=(0,2))
                data_3d = 0.25 * (rot_1 + rot_2 + rot_3 + rot_4)
            elif plane == "z":
                rot_1 = np.rot90(data_3d, axes=(0,1))
                rot_2 = np.rot90(rot_1  , axes=(0,1))
                rot_3 = np.rot90(rot_2  , axes=(0,1))
                rot_4 = np.rot90(rot_3  , axes=(0,1))
                data_3d = 0.25 * (rot_1 + rot_2 + rot_3 + rot_4)
        return data_3d.flatten()


    def filter(self, values):
        # Reshape the linearized data into a 3D array
        data_3d = np.reshape(values, (self.width, self.height, self.length))

        # Perform your filtering operation here
        # Example: Apply a simple 3D Gaussian filter
        sigma = self.proximity_radius  # Adjust the value as needed
        filtered_data = gaussian_filter(data_3d, sigma=sigma)

        # Return the filtered 3D array
        return filtered_data.flatten()

    def plot_cuboid(self, densities, threshold=0.5):
        plot.plot_cuboid(self, densities=densities, threshold=threshold)

    def set_support(self, id, x=None, y=None, z=None):
        self.supports[id] = (x,y,z)

    def set_exponent(self, exponent):
        self.exponent = exponent

    def set_force(self, id, x=0, y=0, z=0):
        self.loads[id] = (x,y,z)

    def set_solver(self, use_cpu=True, use_direct=True):
        self.cpu = use_cpu
        self.direct = use_direct

    def set_path(self, path):
        self.path = path

    def write(self, densities, file):
        with open(file, "w") as out:
            out.write(write.generate_model(self, densities))

    def write_stl(self, file, densities, threshold=0.5):
        write.mesh_to_stl(self.nodes, self.elements, densities, file, threshold=threshold)

    def write_ply(self, file, densities, threshold=0.5):
        write.mesh_to_ply(self.nodes, self.elements, densities, file, threshold=threshold)

    def run(self, densities):
        return run.run(self, densities)

    def _generate_model(self):
        return geom.generate_geometry(width=self.width, height=self.height, length=self.length)

    def _compute_element_midpoints(self):
        return geom.compute_element_midpoints(self.nodes, self.elements)

    def save_model(self, filename):
        with open(filename, 'wb') as file:
            pickle.dump(self, file)

    @staticmethod
    def load_model(filename):
        with open(filename, 'rb') as file:
            model = pickle.load(file)
        return model