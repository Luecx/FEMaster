
class ThermalExpansionIsotropic:
    def __init__(self, alpha: float):
        self.alpha = float(alpha)

    def to_femaster(self):
        return ["*EXPANSION", f"{self.alpha}"]
