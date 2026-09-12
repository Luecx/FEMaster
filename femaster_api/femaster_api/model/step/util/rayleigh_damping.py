"""Mass- and stiffness-proportional Rayleigh damping coefficients.

Rayleigh damping is a transient-step control.  The object stores the two
coefficients explicitly so the surrounding step can emit the native damping
block without a generic options dictionary.
"""

from __future__ import annotations


class RayleighDamping:
    """Rayleigh damping coefficients ``alpha`` and ``beta``."""

    def __init__(self, alpha: float = 0.0, beta: float = 0.0) -> None:
        self.alpha = float(alpha)
        self.beta = float(beta)
