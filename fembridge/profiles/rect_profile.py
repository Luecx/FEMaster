from __future__ import annotations
import math
from .profile_base import Profile


class RectProfile(Profile):
    """
    Rectangle of width b (along y) and height h (along z), centered at origin.
    Centroid at (0,0). Shear center coincides with centroid for a rectangle.

    Parameters
    ----------
    name : str
        Name of the profile (used as reference in *BEAM SECTION).
    b : float
        Width along local y-axis [m].
    h : float
        Height along local z-axis [m].
    """

    def __init__(self, name: str, b: float, h: float):
        b = float(b)
        h = float(h)

        A  = b * h
        Iy = b * h**3 / 12.0  # bending about local y-axis (z^2 distribution)
        Iz = h * b**3 / 12.0  # bending about local z-axis (y^2 distribution)

        # torsion constant (Saint-Venant) — approximate engineering formula
        # It ≈ b*h^3*(1/3 - 0.21*b/h*(1 - b^4/(12*h^4))) for h >= b
        # symmetrized for arbitrary b/h
        if h >= b:
            It = b * h**3 * (1/3 - 0.21 * (b / h) * (1 - (b**4) / (12 * h**4)))
        else:
            It = h * b**3 * (1/3 - 0.21 * (h / b) * (1 - (h**4) / (12 * b**4)))

        ey = 0.0  # shear center offset y
        ez = 0.0  # shear center offset z

        super().__init__(name, Iy, Iz, It, A, ey, ez)
