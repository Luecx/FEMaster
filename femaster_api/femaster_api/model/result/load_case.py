"""One result load case."""

from ..field.field import Field
from .frame import Frame


class LoadCase:
    """Ordered collection of result frames for one load case."""

    def __init__(self, id: int, name: str | None = None) -> None:
        self.id = int(id)
        self.name = name
        self.frames: list[Frame] = []

    def add(self, frame: Frame) -> Frame:
        if any(item.id == frame.id for item in self.frames):
            raise ValueError(f"duplicate frame id: {frame.id}")
        self.frames.append(frame)
        return frame

    def frame(self, key: int | str = 0) -> Frame:
        for frame in self.frames:
            if isinstance(key, int) and frame.id == key:
                return frame
            if isinstance(key, str) and frame.name == key:
                return frame
        raise KeyError(f"unknown frame: {key}")

    def field(self, name: str, *, frame: int | str = 0) -> Field:
        return self.frame(frame).field(name)
