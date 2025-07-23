import numpy as np

class Plane:
    def __init__(self, point, normal):
        """
        Initialize a plane from a base point and a normal vector.

        Parameters
        ----------
        point : array-like, shape (3,)
            A point on the plane.
        normal : array-like, shape (3,)
            A vector normal to the plane.
        """
        self.point = np.asarray(point, dtype=np.float64)
        self.normal = np.asarray(normal, dtype=np.float64)
        self.normal /= np.linalg.norm(self.normal)

    def signed_distance(self, points):
        """
        Compute signed distances from one or more points to the plane.

        Parameters
        ----------
        points : ndarray, shape (..., 3)
            Points to compute the signed distance from.

        Returns
        -------
        distances : ndarray
            Signed distances from the plane (positive above, negative below).
        """
        points = np.atleast_2d(points)
        return (points - self.point) @ self.normal

    def orthonormal_basis(self):
        """
        Construct two orthonormal vectors (r, s) lying in the plane.

        Returns
        -------
        r : ndarray, shape (3,)
            First orthonormal vector in the plane.
        s : ndarray, shape (3,)
            Second orthonormal vector in the plane, perpendicular to r and to the normal.
        """
        # Find any vector not colinear with the normal
        if abs(self.normal[0]) < 0.9:
            v = np.array([1.0, 0.0, 0.0])
        else:
            v = np.array([0.0, 1.0, 0.0])

        r = np.cross(self.normal, v)
        r /= np.linalg.norm(r)
        s = np.cross(self.normal, r)
        return r, s

    def all_on_one_side(self, points, tol=1e-6):
        """
        Check if all points lie on the same side or within tolerance of the plane.

        Parameters
        ----------
        points : ndarray, shape (..., 3)
            Points to check.
        tol : float
            Tolerance for being "on" the plane.

        Returns
        -------
        bool
        """
        dists = self.signed_distance(points)
        return np.all(dists >= -tol) or np.all(dists <= tol)
