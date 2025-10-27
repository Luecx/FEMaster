# sections/point_mass_section.py
from __future__ import annotations
from typing import Tuple
from .section_base import Section
from ..sets.nodeset import NodeSet


class PointMassSection(Section):
    """
    Point mass + (optional) lumped inertias and springs at a NodeSet.

    FEMaster serialization (4 lines after header):
      *POINT MASS, NSET=<name>
      MASS
      INERTIA_X, INERTIA_Y, INERTIA_Z
      Kx, Ky, Kz
      Krx, Kry, Krz
    """
    def __init__(
            self,
            nset       : NodeSet,
            mass       : float,
            inertias   : Tuple[float, float, float] = (0.0, 0.0, 0.0),
            springs    : Tuple[float, float, float] = (0.0, 0.0, 0.0),
            rot_springs: Tuple[float, float, float] = (0.0, 0.0, 0.0),
    ):
        self.nset        = nset
        self.mass        = float(mass)
        self.inertias    = tuple(float(v) for v in inertias)
        self.springs     = tuple(float(v) for v in springs)
        self.rot_springs = tuple(float(v) for v in rot_springs)

        if len(self.inertias) != 3 or len(self.springs) != 3 or len(self.rot_springs) != 3:
            raise ValueError("inertias, springs, and rot_springs must be 3-tuples")

    def to_femaster(self) -> str:
        ix, iy, iz    = self.inertias
        kx, ky, kz    = self.springs
        krx, kry, krz = self.rot_springs
        return "\n".join([
            f"*POINT MASS SECTION, NSET={self.nset.name}",
            f"{self.mass}",
            f"{ix}, {iy}, {iz}",
            f"{kx}, {ky}, {kz}",
            f"{krx}, {kry}, {krz}",
        ])

    def to_asami(self) -> str:
        raise NotImplementedError
