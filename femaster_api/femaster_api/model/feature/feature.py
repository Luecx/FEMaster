"""Base non-topological model feature."""


class Feature:
    """Base class for one assembly-level non-topological feature."""

    def export(self) -> str:
        raise NotImplementedError
