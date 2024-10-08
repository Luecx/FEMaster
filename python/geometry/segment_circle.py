"""
################################################################################
# circle_segment.py
#
# Implementation of a circular arc segment.
#
# Inherits from Segment and uses start and end points with a direction.
#
# Author: Finn Eggers
# Date: 2024-10-08
################################################################################
"""

import numpy as np
from .segment import Segment

class CircleSegment(Segment):
    def __init__(self, start, end, direction, subdivisions, name=None):
        """
        Initialize a circular arc segment.

        Parameters:
        -----------
        direction : int
            Direction of the arc (1 for counterclockwise, -1 for clockwise).
        """
        super().__init__(start, end, subdivisions, name)
        self.direction = direction

    def get_points(self):
        """
        Generates a circular arc segment.

        Returns:
        --------
        points : list of [float, float]
            List of coordinates defining the arc segment.
        """
        # Compute center and radius
        mid_x = (self.start[0] + self.end[0]) / 2
        mid_y = (self.start[1] + self.end[1]) / 2
        radius = np.sqrt((self.end[0] - mid_x)**2 + (self.end[1] - mid_y)**2)

        # Angle between start and end
        start_angle = np.arctan2(self.start[1] - mid_y, self.start[0] - mid_x)
        end_angle = np.arctan2(self.end[1] - mid_y, self.end[0] - mid_x)

        if self.direction == -1:
            start_angle, end_angle = end_angle, start_angle

        # Generate points along the arc
        angles = np.linspace(start_angle, end_angle, self.subdivisions + 1)
        points = [[mid_x + radius * np.cos(angle), mid_y + radius * np.sin(angle)] for angle in angles]
        return points
