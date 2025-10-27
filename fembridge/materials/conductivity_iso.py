
class ConductivityIsotropic:
    def __init__(self, k: float):
        self.k = float(k)

    def to_femaster(self):
        return ["*CONDUCTIVITY", f"{self.k}"]
