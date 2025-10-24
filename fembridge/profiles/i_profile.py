# profiles/i_profile.py
from __future__ import annotations
from .profile_base import Profile


class IProfile(Profile):
    """
    Symmetric I/H-section with:
      h  : overall height (z-direction)
      bf : flange width (y-direction, same top/bottom)
      tf : flange thickness
      tw : web thickness (web height = h - 2*tf)

    Local axes follow your convention (width along y, height along z):
      Iy -> bending about y (z^2 distribution)  ~ O(h^3)
      Iz -> bending about z (y^2 distribution)  ~ O(bf^3)

    Uses standard composite-rectangle formulas; centroid at origin; shear center = centroid.
    Torsion constant It uses a common engineering approximation for thin-walled open sections.
    """
    def __init__(self, name: str, h: float, bf: float, tf: float, tw: float):
        h  = float(h); bf = float(bf); tf = float(tf); tw = float(tw)

        if tf <= 0.0 or tw <= 0.0 or bf <= 0.0 or h <= 0.0:
            raise ValueError("All dimensions must be positive (h, bf, tf, tw > 0).")

        hw = h - 2.0 * tf  # clear web height
        if hw <= 0.0:
            raise ValueError("Web height must be positive: require h > 2*tf.")

        # Area
        A = 2.0 * bf * tf + tw * hw

        # --- Second moments about centroid ---
        # About y (z^2 distribution): strong axis in this convention
        # Top+bottom flange: each has centroidal Iy = bf*tf^3/12, shifted by dz
        dz = (h / 2.0) - (tf / 2.0)  # distance from section centroid to flange centroid (along z)
        Iy_flange_cent = bf * tf**3 / 12.0
        Iy_flange = Iy_flange_cent + bf * tf * dz**2
        Iy_web    = tw * hw**3 / 12.0
        Iy = 2.0 * Iy_flange + Iy_web

        # About z (y^2 distribution)
        Iz_flange = bf**3 * tf / 12.0
        Iz_web    = tw**3 * hw / 12.0
        Iz = 2.0 * Iz_flange + Iz_web

        # --- Saint-Venant torsion constant (engineering approximation) ---
        # J ≈ 2*(bf*tf^3)/3 + (hw*tw^3)/3
        It = 2.0 * (bf * tf**3) / 3.0 + (hw * tw**3) / 3.0

        # Shear center at centroid for doubly symmetric I
        ey = 0.0
        ez = 0.0

        super().__init__(name, Iy, Iz, It, A, ey, ez)
