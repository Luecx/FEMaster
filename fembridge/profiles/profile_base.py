from __future__ import annotations


class Profile:
    """
    Base class for geometric beam profiles.

    Each profile defines:
      Iy, Iz, It, A  – geometric properties
      ey, ez          – offset of shear center to center of gravity
    """

    def __init__(self, name: str, Iy: float, Iz: float, It: float, A: float,
                 ey: float = 0.0, ez: float = 0.0, Iyz: float = 0.0):
        self.name = str(name)
        self.Iy  = float(Iy)
        self.Iz  = float(Iz)
        self.Iyz = float(Iyz)
        self.It  = float(It)
        self.A   = float(A)
        self.ey  = float(ey)
        self.ez  = float(ez)

    def to_femaster(self) -> str:
        """
        Serialize as full FEMaster block.
        FEMaster currently does not support shear-center offsets (ey, ez);
        raises an error if they are nonzero and omits them otherwise.
        """
        if abs(self.ey) > 1e-12 or abs(self.ez) > 1e-12:
            raise ValueError(
                f"FEMaster export for profile '{self.name}' does not support "
                f"nonzero shear-center offsets (ey={self.ey}, ez={self.ez})."
            )
        if abs(self.Iyz) > 1e-12:
            raise ValueError(
                f"FEMaster export for profile '{self.name}' does not support "
                f"nonzero product of inertia Iyz={self.Iyz}."
            )

        return (
            f"*PROFILE, NAME={self.name}\n"
            f"{self.Iy}, {self.Iz}, {self.It}, {self.A}"
        )

    def to_asami(self) -> str:
        raise NotImplementedError("Profile.to_asami not implemented yet.")
