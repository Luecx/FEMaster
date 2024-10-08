"""
################################################################################
# segment_group.py
#
# Implementation of the SegmentGroup class, which groups multiple segments.
#
# Inherits from Segment, allowing SegmentGroup to behave like a single segment
# when necessary. It can check if it forms a closed loop and generate a complete
# list of points for all contained segments.
#
# Author: Finn Eggers
# Date: 2024-10-08
################################################################################
"""

from .segment import Segment
import numpy as np

class SegmentGroup(Segment):
    def __init__(self, segments, name=None, tolerance=1e-6):
        """
        Initialize a SegmentGroup with a list of segments.

        Parameters:
        -----------
        segments : list of Segment
            A list of segment objects (e.g., StraightSegment, CurvedSegment, CircleSegment).
        name : str, optional
            Name of the segment group. If None, a unique name will be assigned.
        tolerance : float, optional
            Tolerance for comparing points to determine continuity and closure.
        """
        # Start and end are determined by the first and last segments in the group
        start = segments[0].start if segments else [0, 0]
        end = segments[-1].end if segments else [0, 0]
        subdivisions = sum([seg.subdivisions for seg in segments])

        super().__init__(start, end, subdivisions, name)
        self.segments = segments
        self.closed = self.is_closed(tolerance)
        self.tolerance = tolerance

    def is_closed(self, tolerance=1e-6):
        """
        Determine if the segment group forms a closed loop.

        Parameters:
        -----------
        tolerance : float
            Tolerance for comparing points to determine if the first and last points are the same.

        Returns:
        --------
        closed : bool
            True if the first segment's start point is within tolerance of the last segment's end point.
        """
        if not self.segments:
            return False
        return self._points_within_tolerance(self.segments[0].start, self.segments[-1].end, tolerance)

    def _points_within_tolerance(self, point1, point2, tolerance):
        """
        Check if two points are within a given tolerance.

        Parameters:
        -----------
        point1 : list or tuple of float
            The first point [x, y].
        point2 : list or tuple of float
            The second point [x, y].
        tolerance : float
            The tolerance value for comparison.

        Returns:
        --------
        within_tolerance : bool
            True if the points are within the specified tolerance, False otherwise.
        """
        return np.linalg.norm(np.array(point1) - np.array(point2)) <= tolerance

    def get_points(self):
        """
        Generate a complete list of points for the entire segment group.

        Returns:
        --------
        points : list of [float, float]
            List of coordinates defining all segments in order.
        """
        points = []
        for segment in self.segments:
            segment_points = segment.get_points()
            # Check if the current segment start is close enough to the last added point
            if points and self._points_within_tolerance(segment_points[0], points[-1], self.tolerance):
                points.extend(segment_points[1:])
            else:
                points.extend(segment_points)

        if self.closed:
            points[-1] = points[0]

        return points
