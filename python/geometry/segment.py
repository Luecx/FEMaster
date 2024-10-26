"""
################################################################################
# segment.py
#
# Abstract base class for different segment types in boundary geometry definition.
#
# Each segment has a start and end point, a number of subdivisions, and an optional name.
# If a name is not provided, it will be automatically assigned a unique name.
#
# Author: Finn Eggers
# Date: 2024-10-08
################################################################################
"""

import numpy as np

class Segment:
    _counter = 0

    def __init__(self, function, t_start=0, t_end=1, n_subdivisions=None, n_points=None, name=None):

        # function to evaluate the segment
        self.func = function

        # starting points and t-values coreesponding to the start and end points
        self.t_start = None
        self.t_end = None
        self.p_start = None
        self.p_end = None

        # set the variables above
        self.truncate(t_start, t_end)

        # number of subdivisions and points
        self.set_detail(n_subdivisions, n_points)

        # name and id
        Segment._counter += 1
        self.id = Segment._counter

        # Assign a unique name if none is provided
        if name is None:
            self.name = f"Segment_{self.id }"
        else:
            self.name = name

    @staticmethod
    def equal_points(point1, point2, tolerance=1e-6):
        return np.linalg.norm(np.array(point1) - np.array(point2)) <= tolerance

    @staticmethod
    def inbetween_points(point, point1, point2, tolerance=1e-6):
        # Convert points to numpy arrays
        p = np.array(point)
        p1 = np.array(point1)
        p2 = np.array(point2)

        # Calculate the vectors
        v10 = p - p1  # Vector from point1 to point
        v12 = p2 - p1  # Vector from point1 to point2

        # Check if the point is along the line segment within tolerance
        proj_length = np.dot(v10, v12) / np.dot(v12, v12)  # Projection of v10 on v12
        if proj_length < -tolerance or proj_length > 1 + tolerance:
            return False  # Not in the segment

        # Check if the point is within the tolerance distance from the line segment
        closest_point = p1 + proj_length * v12
        distance = np.linalg.norm(p - closest_point)

        return distance <= tolerance

    def contains(self, point, tolerance=1e-6):
        points = self.get_points()
        for i in range(len(points) - 1):
            if self.inbetween_points(point, points[i], points[i+1], tolerance):
                return True
        return False

    def get_points(self):
        points = [self.at(t) for t in np.linspace(self.t_start, self.t_end, self.n_points)]
        points[0]  = self.p_start
        points[-1] = self.p_end
        return points

    def at(self, t):
        return self.func(t)

    def truncate(self, t_start=0, t_end=1):
        self.t_start = t_start
        self.t_end = t_end

        self.p_start = np.asarray(self.at(t_start))
        self.p_end   = np.asarray(self.at(t_end))

    def set_function(self, function):
        self.func = function
        self.truncate(self.t_start, self.t_end)

    def set_detail(self, n_subdivisions=None, n_points=None):
        if n_subdivisions is not None and n_points is not None:
            raise ValueError("Only one of 'n_subdivisions' or 'n_points' can be provided.")
        if n_subdivisions is None and n_points is None:
            raise ValueError("Either 'n_subdivisions' or 'n_points' must be provided.")
        if n_subdivisions is not None:
            self.n_subdivisions = n_subdivisions
            self.n_points = n_subdivisions + 1
        else:
            self.n_points = n_points
            self.n_subdivisions = n_points - 1

        if self.n_subdivisions < 1:
            self.n_subdivisions = 1
            self.n_points = 2

    def tangent(self, t):
        dir = self.at(t + 1e-6) - self.at(t - 1e-6)
        return dir / np.linalg.norm(dir)

    def perpendicular(self, t):
        tan = self.tangent(t)
        return np.array([-tan[1], tan[0]])

    def plot(self, ax=None, **kwargs):
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()

        points = self.get_points()
        x, y = zip(*points)
        ax.plot(x, y, **kwargs)
        ax.scatter(x, y, color='black')
        plt.show()
