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
    def __init__(self, start, end, n_subdivisions=None, name=None):
        """
        Initialize a straight line segment with start and end points.

        Parameters:
        -----------
        start : list or tuple of float
            Coordinates of the start point [x, y].
        end : list or tuple of float
            Coordinates of the end point [x, y].
        """
        # Define the linear function for the segment
        def linear_function(t):
            return np.array([
                start[0] + t * (end[0] - start[0]),
                start[1] + t * (end[1] - start[1])
            ])

        # Initialize the base class with the linear function
        super().__init__(linear_function, t_start=0, t_end=1, n_subdivisions=n_subdivisions, name=name)

        # Store the start and end points for reference
        self.start = start
        self.end = end
