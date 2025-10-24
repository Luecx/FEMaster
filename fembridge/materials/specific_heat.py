
class SpecificHeat:
    def __init__(self, cp: float):
        self.cp = float(cp)

    def to_femaster_lines(self):
        return ["*SPECIFIC HEAT", f"{self.cp}"]
