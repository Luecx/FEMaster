"""Base assembly-level constraint."""


class Constraint:
    """Base class for one assembly-level constraint."""

    def export(self) -> str:
        raise NotImplementedError
