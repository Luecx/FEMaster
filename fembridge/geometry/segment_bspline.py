import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
from .segment import Segment

class BSplineSegment(Segment):
    def __init__(self, control_points, closed=False, n_subdivisions=1000, name=None):
        """
        Initialize a B-spline segment with control points.

        Parameters:
        -----------
        control_points : list or ndarray of shape (N, 2)
            List of 2D control points defining the B-spline.
        closed : bool, optional
            If True, forces the B-spline to be closed. Default is False.
        n_subdivisions : int, optional
            Number of subdivisions for discretizing the curve. Default is 1000.
        name : str, optional
            Name of the segment.
        """
        self.control_points = np.array(control_points)
        self.closed = closed
        self.n_subdivisions = n_subdivisions
        self.name = name

        if self.closed:
            # Append the starting coordinates to close the curve
            x = np.r_[self.control_points[:, 0], self.control_points[0, 0]]
            y = np.r_[self.control_points[:, 1], self.control_points[0, 1]]
            # Fit splines to x=f(u) and y=g(u), treating both as periodic
            self.tck, self.u = interpolate.splprep([x, y], s=0, per=True)
        else:
            x = self.control_points[:, 0]
            y = self.control_points[:, 1]
            # Fit splines to x=f(u) and y=g(u) without periodic constraint
            self.tck, self.u = interpolate.splprep([x, y], s=0, per=False)

        # Define the B-spline function for evaluation
        def bspline_function(t):
            xi, yi = interpolate.splev(t, self.tck)
            return np.array([xi, yi]).T

        # Initialize the base Segment class
        super().__init__(bspline_function, t_start=0, t_end=1, n_subdivisions=n_subdivisions, name=name)

    def is_closed(self):
        """Return True if the B-spline is closed."""
        return self.closed