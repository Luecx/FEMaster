"""
################################################################################
# straight_segment.py
#
# Implementation of a straight line segment.
#
# Inherits from Segment and generates equally spaced points between start and end.
#
# Author: Finn Eggers
# Date: 2024-10-08
################################################################################
"""

import numpy as np
from .segment import Segment

class StraightSegment(Segment):
    def get_points(self):
        """
        Generates equally spaced points between the start and end points.

        Returns:
        --------
        points : list of [float, float]
            List of coordinates along the straight segment.
        """
        x_vals = np.linspace(self.start[0], self.end[0], self.subdivisions + 1)
        y_vals = np.linspace(self.start[1], self.end[1], self.subdivisions + 1)
        return [[x, y] for x, y in zip(x_vals, y_vals)]
