import numpy as np
import geom
import plot
import filter
import write
import subprocess
import os
import run

class Model:
    def __init__(self, width, height, length):
        self.width  = width
        self.height = height
        self.length = length
        # geometry
        self.nodes, self.elements = self._generate_model()
        self.midpoints = self._compute_element_midpoints()

        # material
        self.youngs = 210000
        self.nu = 0.3

        # filtering
        self.proximities = []
        self.set_proximity_radius(width / 5)

        # constraints
        self.supports = {}
        self.loads    = {}

        # simp exonent
        self.exponent = 3

        # solver mode
        self.cpu = True
        self.direct = True

    def get_special_node_ids(self):
        return geom.get_special_node_ids(width=self.width, height=self.height, length=self.length)

    def set_proximity_radius(self, radius):
        self.proximities = []
        for i, midpoint in enumerate(self.midpoints):
            close_elements = [j for j, m in enumerate(self.midpoints) if np.linalg.norm(np.array(midpoint) - np.array(m)) <= radius and j != i]
            self.proximities.append(close_elements)

    def set_material(self, youngs, nu):
        self.youngs = youngs
        self.nu = nu

    def filter(self, values):
        return filter.filter(values, self.proximities)

    def plot_cuboid(self, densities, threshold=0.5):
        plot.plot_cuboid(self.nodes, self.elements, densities=densities, threshold=threshold)

    def set_support(self, id, x=None, y=None, z=None):
        self.supports[id] = (x,y,z)

    def set_exponent(self, exponent):
        self.exponent = exponent

    def set_force(self, id, x=0, y=0, z=0):
        self.loads[id] = (x,y,z)

    def set_solver(self, use_cpu=True, use_direct=True):
        self.cpu = use_cpu
        self.direct = use_direct

    def write(self, densities, file):
        with open(file, "w") as out:
            out.write(write.generate_model(self, densities))

    def run(self, densities):
        return run.run(self, densities)

    def _generate_model(self):
        return geom.generate_geometry(width=self.width, height=self.height, length=self.length)

    def _compute_element_midpoints(self):
        return geom.compute_element_midpoints(self.nodes, self.elements)