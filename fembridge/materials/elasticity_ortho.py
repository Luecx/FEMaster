
class ElasticityOrthotropic:
    def __init__(self, *, Ex: float, Ey: float, Ez: float,
                 nuxy: float, nuxz: float, nyz: float,
                 Gxy: float, Gxz: float, Gyz: float):
        self.Ex = float(Ex)
        self.Ey = float(Ey)
        self.Ez = float(Ez)

        self.nuxy = float(nuxy)
        self.nuxz = float(nuxz)
        self.nyz  = float(nyz)

        self.Gxy = float(Gxy)
        self.Gxz = float(Gxz)
        self.Gyz = float(Gyz)

    def to_femaster(self):
        return [
            "*ELASTIC, TYPE=ORTHOTROPIC",
            f"{self.Ex}, {self.Ey}, {self.Ez}, {self.nuxy}, {self.nuxz}, {self.nyz}, {self.Gxy}, {self.Gxz}, {self.Gyz}",
        ]
