"""Format-independent solver result."""

from ..field.field import Field
from .load_case import LoadCase


class Result:
    """Complete format-independent result of one FEMaster execution."""

    def __init__(self) -> None:
        self.loadcases: list[LoadCase] = []

    def add(self, loadcase: LoadCase) -> LoadCase:
        if any(item.id == loadcase.id for item in self.loadcases):
            raise ValueError(f"duplicate loadcase id: {loadcase.id}")
        self.loadcases.append(loadcase)
        return loadcase

    def loadcase(self, key: int | str = 1) -> LoadCase:
        for loadcase in self.loadcases:
            if isinstance(key, int) and loadcase.id == key:
                return loadcase
            if isinstance(key, str) and loadcase.name == key:
                return loadcase
        raise KeyError(f"unknown loadcase: {key}")

    def field(
        self,
        name: str,
        *,
        loadcase: int | str = 1,
        frame: int | str = 0,
    ) -> Field:
        return self.loadcase(loadcase).field(name, frame=frame)
