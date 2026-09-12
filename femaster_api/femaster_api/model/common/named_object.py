"""Immutable semantic-name base object."""


class NamedObject:
    """Base object with an immutable non-empty semantic name."""

    __slots__ = ("_name",)

    def __init__(self, name: str) -> None:
        name = str(name).strip()
        if not name:
            raise ValueError("name must not be empty")
        self._name = name

    @property
    def name(self) -> str:
        return self._name

    def __repr__(self) -> str:
        return f"{type(self).__name__}(name={self.name!r})"
