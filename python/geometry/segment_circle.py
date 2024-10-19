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
    def __init__(self, center, start_point, angles, subdivisions, name=None):
        """
        Initialize a circular arc segment.

        Parameters:
        -----------
        direction : int
            Direction of the arc (1 for counterclockwise, -1 for clockwise).
        """
        self.center = center
        self.start_angle = np.arctan2(start_point[1] - center[1], start_point[0] - center[0])
        self.end_angle   = self.start_angle + angles
        self.angles      = angles
        self.radius      = np.sqrt((start_point[0] - self.center[0]) ** 2 + (start_point[1] - self.center[1]) ** 2)
        self.end_point   = [self.center[0] + np.cos(self.end_angle) * self.radius,
                            self.center[1] + np.sin(self.end_angle) * self.radius]

        super().__init__(start_point, self.end_point, subdivisions, name)

    def get_points(self):
        """
        Generates a circular arc segment.

        Returns:
        --------
        points : list of [float, float]
            List of coordinates defining the arc segment.
        """
        angles = np.linspace(self.start_angle, self.end_angle, self.subdivisions)
        # remove first and last
        angles = angles[1:-1]

        # compute points for each angle
        points = [[self.center[0] + self.radius * np.cos(k),
                   self.center[1] + self.radius * np.sin(k)] for k in angles]

        # append the start and end point
        points.insert(0, self.start)
        points.append(self.end)
        return points
