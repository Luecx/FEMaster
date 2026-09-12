"""Base collector-owned load."""


class Load:
    """Base class for one collector-owned load definition."""

    def export(self, collector: str) -> str:
        raise NotImplementedError
