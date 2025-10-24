from __future__ import annotations

class Section:
    """Base class for all sections. No IDs or names; purely applied data."""
    def to_femaster(self) -> str:
        raise NotImplementedError

    def to_asami(self) -> str:
        raise NotImplementedError
