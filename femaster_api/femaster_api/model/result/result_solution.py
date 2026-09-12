"""One solver/loadcase solution inside a result file.

A ``Solution`` corresponds to one FEMaster loadcase or FRD analysis step and owns
an ordered sequence of output frames.  The optional ``type`` stores the procedure
category when the source result format exposes it, for example ``STATIC``,
``DYNAMIC``, ``MODAL`` or ``BUCKLING``.

Frame IDs are unique within a solution.  Physical time/frequency/load values live
on the frames themselves rather than being encoded in field names.
"""

from __future__ import annotations

from ..field.field import Field
from ..field.field_type import FieldType
from .result_frame import Frame


class Solution:
    """Ordered collection of frames produced by one analysis solution."""

    def __init__(
        self,
        id: int,
        *,
        name: str | None = None,
        type: str | None = None,
    ) -> None:
        self.id = int(id)
        self.name = name
        self.type = type
        self.frames: list[Frame] = []

    def add(self, frame: Frame) -> Frame:
        """Append one frame while enforcing a unique frame ID."""

        if any(item.id == frame.id for item in self.frames):
            raise ValueError(
                f"duplicate frame id {frame.id} in solution {self.id}"
            )
        self.frames.append(frame)
        return frame

    def frame(self, key: int | str = 0) -> Frame:
        """Return a frame by numeric ID or optional semantic name."""

        for frame in self.frames:
            if isinstance(key, int) and frame.id == key:
                return frame
            if isinstance(key, str) and frame.name == key:
                return frame
        raise KeyError(f"unknown frame in solution {self.id}: {key}")

    def field(
        self,
        key: str | FieldType,
        *,
        frame: int | str = 0,
    ) -> Field:
        """Return one field from a selected frame."""

        return self.frame(frame).field(key)

    def __len__(self) -> int:
        return len(self.frames)
