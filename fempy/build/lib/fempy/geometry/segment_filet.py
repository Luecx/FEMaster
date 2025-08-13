from .segment import Segment
from scipy.optimize import minimize
import numpy as np

class Filet(Segment):
    def __init__(self, segment1, segment2, radius, n_subdivisions=None, name=None):
        self.segment1 = segment1
        self.segment2 = segment2
        self.radius = radius

        # Get optimized t-values on the two segments
        t1, t2 = self._optimize_t_values()
        self.t1 = t1
        self.t2 = t2

        # Calculate points and tangent directions at optimized t-values
        p1 = segment1.at(t1)
        p2 = segment2.at(t2)
        tangent1 = segment1.tangent(t1)
        tangent2 = segment2.tangent(t2)

        # Determine the center of the circular arc
        arc_center = self._calculate_arc_center(p1, p2, tangent1, tangent2)

        # Define the circular arc function as the filet shape
        filet_function = lambda t: self._arc_point(arc_center, p1, p2, t)

        # Initialize Segment with the arc as the filet
        super().__init__(filet_function, 0, 1, n_subdivisions=n_subdivisions, name=name)

    def _optimize_t_values(self):
        def objective(t_values):
            # Compute points and tangents
            point1 = self.segment1.at(t_values[0])
            point2 = self.segment2.at(t_values[1])
            tangent1 = self.segment1.tangent(t_values[0])
            tangent2 = self.segment2.tangent(t_values[1])

            # Compute arc center
            arc_center = self._calculate_arc_center(point1, point2, tangent1, tangent2)

            # Minimize distance between arc center and both points to get best fit
            dist1 = np.linalg.norm(arc_center - point1) - self.radius
            dist2 = np.linalg.norm(arc_center - point2) - self.radius

            # Return sum of squared residuals
            return dist1 ** 2 + dist2 ** 2

        def newton_objective(t_values):
            # Compute residuals for Newton method
            point1 = self.segment1.at(t_values[0])
            point2 = self.segment2.at(t_values[1])
            tangent1 = self.segment1.tangent(t_values[0])
            tangent2 = self.segment2.tangent(t_values[1])
            arc_center = self._calculate_arc_center(point1, point2, tangent1, tangent2)
            dist1 = np.linalg.norm(arc_center - point1) - self.radius
            dist2 = np.linalg.norm(arc_center - point2) - self.radius
            return np.array([dist1, dist2])

        import random
        from scipy.optimize import minimize, root

        # Initial guesses and bounds
        initial_guess = [self.segment1.t_end, self.segment2.t_start]
        bounds = [(self.segment1.t_start, self.segment1.t_end),
                  (self.segment2.t_start, self.segment2.t_end)]
        initial_guesses = [initial_guess] + [
            [random.uniform(bounds[0][0], bounds[0][1]),
             random.uniform(bounds[1][0], bounds[1][1])]
            for _ in range(1, 10)
        ]

        for g in initial_guesses:
            result = minimize(objective, g, bounds=bounds, method='SLSQP')
            if result.success and result.fun < 1e-6:
                t_values = result.x  # Optimized t-values from initial optimization

                # Refine with Newton-Raphson method using `scipy.optimize.root`
                newton_result = root(newton_objective, t_values, method='hybr')

                if newton_result.success:
                    return newton_result.x  # Converged t-values for both segments

        raise ValueError("Filet optimization failed to converge.")

    def _calculate_arc_center(self, p1, p2, tangent1, tangent2):
        # Perpendicular directions to tangents at p1 and p2
        perp1 = np.array([-tangent1[1], tangent1[0]])
        perp2 = np.array([-tangent2[1], tangent2[0]])

        # System of equations to solve for t and s in line equations
        A = np.array([perp1, -perp2]).T
        b = p2 - p1
        t_s = np.linalg.solve(A, b)  # Direct solve instead of least squares

        # Calculate the arc center at the intersection point
        arc_center = p1 + t_s[0] * perp1
        return arc_center


    def _arc_point(self, center, p1, p2, t):
        # Calculate angles for p1 and p2 relative to the center
        angle_start = np.arctan2(p1[1] - center[1], p1[0] - center[0])
        angle_end   = np.arctan2(p2[1] - center[1], p2[0] - center[0])

        if angle_end < angle_start:
            angle_end += 2 * np.pi  # Ensure counterclockwise sweep

        # Interpolate angle for point along the arc
        angle = angle_start + t * (angle_end - angle_start)
        return np.array([center[0] + self.radius * np.cos(angle), center[1] + self.radius * np.sin(angle)])
