import numpy as np

class Mesher:
    def __init__(self, init_size, min_size, curvature_factor):
        self.init_size = init_size
        self.min_size = min_size
        self.curvature_factor = curvature_factor

    def compute(self, geometry):
        length    = geometry.length()
        max_curvature = geometry.curvature_max()

        if max_curvature == 0:
            return self.init_size

        max_curvature_radius = 1 / max_curvature

