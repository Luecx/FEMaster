"""Newmark time-integration parameters for transient analysis.

The two parameters are kept in a dedicated utility object because they describe
the integration algorithm rather than the structural model itself.
"""

from __future__ import annotations


class NewmarkControl:
    """Newmark ``beta`` and ``gamma`` integration parameters."""

    def __init__(self, beta: float = 0.25, gamma: float = 0.5) -> None:
        self.beta = float(beta)
        self.gamma = float(gamma)
