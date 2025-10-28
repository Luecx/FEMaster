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

    def __init__(self, function, t_start=0, t_end=1, n_subdivisions=1, name=None):

        # function to evaluate the segment
        self.func = function

        # starting points and t-values coreesponding to the start and end points
        self.t_start = None
        self.t_end = None
        self.p_start = None
        self.p_end = None

        # number of subdivisions and points
        self.n_subdivisions = None
        self.n_points = None

        # set the variables above
        self.truncate(t_start, t_end)

        # number of subdivisions and points
        self.set_detail(n_subdivisions)

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
        if self.n_subdivisions is None:
            raise ValueError("Number of subdivisions not set.")

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

    def set_detail(self, n_subdivisions=None, mesh_size=None):
        self.n_subdivisions = n_subdivisions
        if n_subdivisions is not None:
            self.n_points = n_subdivisions+1
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

    def length(self):
        length = 0
        x_vals = np.linspace(self.t_start, self.t_end, 1000)

        for i in range(1, len(x_vals)):
            x = self.at(x_vals[i])
            x_prev = self.at(x_vals[i-1])
            length += np.linalg.norm(x - x_prev)

        return length

    def curvature(self, t):
        # Define a small interval around the point `t` for sampling
        delta_t = 1e-3 * (self.t_end - self.t_start)  # Small step based on the total range

        # Get points just before, at, and after `t`
        t_prev = max(self.t_start, t - delta_t)
        t_next = min(self.t_end, t + delta_t)

        x_prev = self.at(t_prev)
        x = self.at(t)
        x_next = self.at(t_next)

        # Calculate vectors between the points
        v1 = x - x_prev
        v2 = x_next - x

        # Calculate the curvature at this single point `t`
        cross_norm = np.linalg.norm(np.cross(v1, v2))
        denominator = np.linalg.norm(v1) * np.linalg.norm(v2) * np.linalg.norm(v1 + v2)

        # If the denominator is zero, return infinity (indicating a straight line)
        if denominator == 0:
            return float('inf')

        # Calculate the curvature radius as the reciprocal of the curvature
        curvature = 2 * cross_norm / denominator
        return curvature

    def min_curvature(self):
        x_vals = np.linspace(self.t_start, self.t_end, 1000)
        curv_min = float('inf')
        for x in x_vals:
            curv_min = min(curv_min, self.curvature(x))
        return curv_min

    def max_curvature(self):
        x_vals = np.linspace(self.t_start, self.t_end, 1000)
        curv_max = 0
        for x in x_vals:
            curv_max = max(curv_max, self.curvature(x))
        return curv_max