
class ElasticityIsotropic:
    def __init__(self, E: float, nu: float):
        self.E = float(E)
        self.nu = float(nu)

    def to_femaster(self):
        return ["*ELASTIC, TYPE=ISOTROPIC", f"{self.E}, {self.nu}"]
