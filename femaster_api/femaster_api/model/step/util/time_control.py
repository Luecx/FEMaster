"""Time interval and step size for transient integration.

The control is intentionally a tiny explicit value object.  It belongs to
``TransientStep`` and does not introduce a separate repository or export
framework.
"""

from __future__ import annotations


class TimeControl:
    """Start time, end time and nominal integration increment."""

    def __init__(self, start: float, end: float, step: float) -> None:
        self.start = float(start)
        self.end = float(end)
        self.step = float(step)
