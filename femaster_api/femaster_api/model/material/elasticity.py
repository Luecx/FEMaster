"""Base elastic constitutive definition."""


class Elasticity:
    """Base class for exportable elastic material behavior."""

    def export(self) -> str:
        raise NotImplementedError
