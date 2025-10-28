"""
################################################################################
# circle_segment.py
#
# Implementation of a circular arc segment.
#
# Inherits from Segment and uses start and end points with a specified angle.
#
# Author: Finn Eggers
# Date: 2024-10-08
################################################################################
"""

import numpy as np
from .segment import Segment

class CircleSegment(Segment):
    def __init__(self, center, start_point, angles, n_subdivisions=None, name=None):
        """
        Initialize a circular arc segment with a center, start point, and angle.

        Parameters:
        -----------
        center : list or tuple of float
            Coordinates of the circle center [x, y].
        start_point : list or tuple of float
            Coordinates of the start point on the circle [x, y].
        angles : float
            Angle to define the arc in radians (positive for counterclockwise).
        """
        # Calculate the initial start and end angles
        self.center = np.array(center)
        self.start_angle = np.arctan2(start_point[1] - center[1], start_point[0] - center[0])
        self.end_angle = self.start_angle + angles
        self.radius = np.linalg.norm(np.array(start_point) - self.center)

        # Define the circular arc function
        def circular_function(t):
            angle = (1 - t) * self.start_angle + t * self.end_angle
            x = self.center[0] + np.cos(angle) * self.radius
            y = self.center[1] + np.sin(angle) * self.radius
            return np.array([x, y])

        # Initialize the base class with the circular function
        super().__init__(circular_function, t_start=0, t_end=1, n_subdivisions=n_subdivisions, name=name)
