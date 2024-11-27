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
    def __init__(self, start, end, mid, n_subdivisions=None, name=None):
        """
        Initialize a curved segment with a start, end, and mid-point.

        Parameters:
        -----------
        mid : list or tuple of float
            Coordinates of the mid-point [x, y].
        """
        # Define the spline function using the quadratic formula with start, mid, and end points
        def spline_function(t):
            x = (1 - t)**2 * start[0] + 2 * (1 - t) * t * mid[0] + t**2 * end[0]
            y = (1 - t)**2 * start[1] + 2 * (1 - t) * t * mid[1] + t**2 * end[1]
            return np.array([x, y])

        # Initialize the base class with the spline function
        super().__init__(spline_function, t_start=0, t_end=1, n_subdivisions=n_subdivisions,name=name)

        # Store the mid-point for reference
        self.mid = mid
