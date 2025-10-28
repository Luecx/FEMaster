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
from .segment_filet import Filet
import numpy as np

class SegmentGroup(Segment):
    def __init__(self, segments, name=None):
        # Store segments and check if the group forms a closed loop
        self.segments   = segments
        self.closed     = self.is_closed()

        at,total_subdivisions = self._create_function()

        # Initialize the Segment base class
        super().__init__(function=at, t_start=0, t_end=1, n_subdivisions=total_subdivisions, name=name)

    def is_closed(self, tolerance=1e-6):
        """
        Check if the segment group forms a closed loop.
        """
        if not self.segments:
            return False
        return Segment.equal_points(self.segments[0].p_start, self.segments[-1].p_end, tolerance)

    def insert_filet(self, seg_id_1, radius, n_subdivisions=5):
        """
        Insert a filet between two segments in the group.
        """
        # Get the segments
        seg_id_2 = (seg_id_1 + 1) % len(self.segments)
        seg1 = self.segments[seg_id_1]
        seg2 = self.segments[seg_id_2]

        filet = Filet(seg1, seg2, radius, n_subdivisions=n_subdivisions)
        ## update t values
        seg1.truncate(seg1.t_start, filet.t1)
        seg2.truncate(filet.t2    , seg2.t_end)


        ## Insert the filet
        self.segments.insert(seg_id_2, filet)
        func, total_subdivisions = self._create_function()
        self.set_function(func)
        self.set_detail(n_subdivisions=total_subdivisions)

    def _create_function(self):
        # Map each segment to its proportional t_start and t_end within the group
        segment_ranges     = []
        total_subdivisions = sum(seg.n_subdivisions for seg in self.segments)
        cumulative_start   = 0
        for seg in self.segments:
            seg_fraction = seg.n_subdivisions / total_subdivisions
            segment_ranges.append((cumulative_start, cumulative_start + seg_fraction))
            cumulative_start += seg_fraction
        self.segment_ranges = segment_ranges

        # functor
        def at(t):
            for (seg, (t_start, t_end)) in zip(self.segments, segment_ranges):
                if t_start -1e-12 <= t <= t_end + 1e-12:
                    local_t = (t - t_start) / (t_end - t_start)
                    local_t = seg.t_start + local_t * (seg.t_end - seg.t_start)
                    return seg.at(local_t)
            return self.segments[-1].at(1)

        return at, total_subdivisions

    def __str__(self):
        """
        print all segments with their id, name, t_start, t_end and class name
        :return:
        """
        s = ""
        for i, seg in enumerate(self.segments):
            s += f"Segment {i}: {seg.name} ({seg.t_start}, {seg.t_end})\n"
        return s