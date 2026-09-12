"""Newmark integration parameters."""


class NewmarkControl:
    """Newmark beta/gamma integration parameters."""

    def __init__(self, beta: float = 0.25, gamma: float = 0.5) -> None:
        self.beta = float(beta)
        self.gamma = float(gamma)
