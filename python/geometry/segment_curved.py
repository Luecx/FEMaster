"""
################################################################################
# curved_segment.py
#
# Implementation of a curved segment using a quadratic spline.
#
# Inherits from Segment and includes a mid-point to define the curvature.
#
# Author: Finn Eggers
# Date: 2024-10-08
################################################################################
"""

import numpy as np
from .segment import Segment

class CurvedSegment(Segment):
    def __init__(self, start, end, mid, subdivisions, name=None):
        """
        Initialize a curved segment with a start, end, and mid-point.

        Parameters:
        -----------
        mid : list or tuple of float
            Coordinates of the mid-point [x, y].
        """
        super().__init__(start, end, subdivisions, name)
        self.mid = mid

    def get_points(self):
        """
        Generates a curved segment using a quadratic spline.

        Returns:
        --------
        points : list of [float, float]
            List of coordinates defining the curved segment.
        """
        t_vals = np.linspace(0, 1, self.subdivisions + 1)
        points = []
        for t in t_vals:
            x = (1 - t)**2 * self.start[0] + 2 * (1 - t) * t * self.mid[0] + t**2 * self.end[0]
            y = (1 - t)**2 * self.start[1] + 2 * (1 - t) * t * self.mid[1] + t**2 * self.end[1]
            points.append([x, y])
        return points
