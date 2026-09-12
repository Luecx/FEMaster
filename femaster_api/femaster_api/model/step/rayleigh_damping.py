"""Rayleigh damping coefficients."""


class RayleighDamping:
    """Mass- and stiffness-proportional damping coefficients."""

    def __init__(self, alpha: float = 0.0, beta: float = 0.0) -> None:
        self.alpha = float(alpha)
        self.beta = float(beta)
