from __future__ import annotations
import math
from .profile_base import Profile

class CircleProfile(Profile):
    """
    Solid circular profile (diameter d). Centroid = shear center.
    """
    def __init__(self, name: str, d: float):
        d = float(d)
        r = d / 2.0

        A  = math.pi * r**2
        Iy = Iz = (math.pi * d**4) / 64.0
        It = (math.pi * d**4) / 32.0  # polar J for solid circle
        ey = 0.0
        ez = 0.0

        super().__init__(name, Iy, Iz, It, A, ey, ez)
