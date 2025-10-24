from __future__ import annotations
import math
from .profile_base import Profile


class TubeProfile(Profile):
    """
    Circular tube (annulus). Provide outer diameter do and wall thickness t.
    Inner diameter di = do - 2t. Centroid = shear center.
    """
    def __init__(self, name: str, do: float, t: float):
        do = float(do)
        t  = float(t)

        di = do - 2.0 * t
        if di <= 0.0:
            raise ValueError("Invalid wall thickness: inner diameter <= 0.")

        A  = math.pi * (do**2 - di**2) / 4.0
        Iy = Iz = math.pi * (do**4 - di**4) / 64.0
        It = math.pi * (do**4 - di**4) / 32.0  # polar moment
        ey = 0.0
        ez = 0.0

        super().__init__(name, Iy, Iz, It, A, ey, ez)
