"""One result frame."""

from ..field.field import Field


class Frame:
    """One result frame inside a load case."""

    def __init__(self, id: int = 0, name: str | None = None) -> None:
        self.id = int(id)
        self.name = name
        self.fields: dict[str, Field] = {}

    def add(self, field: Field) -> Field:
        self.fields[field.name] = field
        return field

    def field(self, name: str) -> Field:
        return self.fields[name]

    def __getitem__(self, name: str) -> Field:
        return self.field(name)
