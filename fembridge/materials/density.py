
class Density:
    def __init__(self, rho: float):
        self.rho = float(rho)

    def to_femaster(self):
        return ["*DENSITY", f"{self.rho}"]
