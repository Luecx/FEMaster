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

from abc import ABC, abstractmethod

import numpy as np

class Segment(ABC):
    _counter = 0  # Class-level counter to generate unique names

    def __init__(self, start, end, subdivisions, name=None):
        """
        Initialize a base segment with start and end points, subdivisions, and an optional name.

        Parameters:
        -----------
        start : list or tuple of float
            Coordinates of the starting point [x, y].
        end : list or tuple of float
            Coordinates of the end point [x, y].
        subdivisions : int
            Number of subdivisions for discretizing the segment.
        name : str, optional
            Name of the segment, useful for node sets or boundary naming.
            If None, a unique name is assigned automatically.
        """
        self.start = start
        self.end = end
        self.subdivisions = subdivisions

        Segment._counter += 1
        self.id = Segment._counter

        # Assign a unique name if none is provided
        if name is None:
            self.name = f"Segment_{self.id }"
        else:
            self.name = name

    @abstractmethod
    def get_points(self):
        """
        Abstract method to be implemented by subclasses to return a list of points
        representing the discretized segment.

        Returns:
        --------
        points : list of [float, float]
            List of coordinates defining the segment.
        """
        pass

    def contains(self, point, tolerance=1e-6):
        """
        Check if a point is approximately on this segment using a piecewise linear approximation.

        Parameters:
        -----------
        point : list or tuple of float
            Coordinates of the point [x, y].
        tolerance : float, optional
            Tolerance for determining if the point is on the segment.

        Returns:
        --------
        is_on_segment : bool
            True if the point is on the segment, False otherwise.
        """
        point = np.array(point)
        segment_points = self.get_points()

        # Check each consecutive pair of points as a straight line segment
        for i in range(len(segment_points) - 1):
            start = np.array(segment_points[i])
            end = np.array(segment_points[i + 1])

            # Compute the vector from start to end and from start to the point
            line_vec = end - start
            point_vec = point - start

            # Check if the cross product is approximately zero (collinearity check)
            if np.abs(np.cross(line_vec, point_vec)) > tolerance:
                continue

            # Check if the point is within the segment bounds using dot product
            if np.dot(point_vec, line_vec) < 0 or np.dot(point_vec, line_vec) > np.dot(line_vec, line_vec):
                continue

            return True

        return False

    def get_id(self):
        return self.id