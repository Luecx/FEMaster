"""Transient time interval control."""


class TimeControl:
    """Start/end/step definition for transient integration."""

    def __init__(self, start: float, end: float, step: float) -> None:
        self.start = float(start)
        self.end = float(end)
        self.step = float(step)
