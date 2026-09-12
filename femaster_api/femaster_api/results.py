"""Format-independent FEMaster result hierarchy.

Readers populate Result -> LoadCase -> Frame -> Field objects. The hierarchy is
shared by RES and FRD readers so downstream code does not depend on file syntax.
"""

from __future__ import annotations

from .fields import Field


class Frame:
    """One result frame inside a load case."""

    def __init__(self, id: int = 0, name: str | None = None) -> None:
        self.id     = int(id)
        self.name   = name
        self.fields: dict[str, Field] = {}

    def add(self, field: Field) -> Field:
        """Add or replace a field by semantic field name."""

        self.fields[field.name] = field
        return field

    def field(self, name: str) -> Field:
        """Return a field by name."""

        return self.fields[name]

    def __getitem__(self, name: str) -> Field:
        return self.field(name)


class LoadCase:
    """Ordered collection of frames produced by one analysis load case."""

    def __init__(self, id: int, name: str | None = None) -> None:
        self.id     = int(id)
        self.name   = name
        self.frames: list[Frame] = []

    def add(self, frame: Frame) -> Frame:
        """Append a result frame while enforcing a unique frame id."""

        if any(item.id == frame.id for item in self.frames):
            raise ValueError(f"duplicate frame id: {frame.id}")
        self.frames.append(frame)
        return frame

    def frame(self, key: int | str = 0) -> Frame:
        """Return a frame by id or optional semantic name."""

        for frame in self.frames:
            if isinstance(key, int) and frame.id == key:
                return frame
            if isinstance(key, str) and frame.name == key:
                return frame
        raise KeyError(f"unknown frame: {key}")

    def field(self, name: str, *, frame: int | str = 0) -> Field:
        """Return one field from a selected frame."""

        return self.frame(frame).field(name)


class Result:
    """Complete format-independent result of one FEMaster execution."""

    def __init__(self) -> None:
        self.loadcases: list[LoadCase] = []

    def add(self, loadcase: LoadCase) -> LoadCase:
        """Append a load case while enforcing a unique numeric id."""

        if any(item.id == loadcase.id for item in self.loadcases):
            raise ValueError(f"duplicate loadcase id: {loadcase.id}")
        self.loadcases.append(loadcase)
        return loadcase

    def loadcase(self, key: int | str = 1) -> LoadCase:
        """Return a load case by id or optional semantic name."""

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
        """Return a field directly from a selected load case and frame."""

        return self.loadcase(loadcase).field(name, frame=frame)
