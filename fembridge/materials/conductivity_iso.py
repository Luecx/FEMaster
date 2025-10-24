
class ConductivityIsotropic:
    def __init__(self, k: float):
        self.k = float(k)

    def to_femaster_lines(self):
        return ["*CONDUCTIVITY", f"{self.k}"]
