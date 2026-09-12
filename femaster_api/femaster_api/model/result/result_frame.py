"""One physical output frame of a FEMaster solution.

A frame groups all fields evaluated at the same analysis state.  ``id`` is the
discrete frame/mode/increment identifier used by the result format, while
``value`` stores the associated physical scalar such as time, frequency,
buckling factor or load parameter when the source format provides it.

Fields are owned directly by the frame and keyed by their native semantic name.
This makes the result hierarchy explicit: ``Result -> Solution -> Frame -> Field``.
"""

from __future__ import annotations

from ..field.field import Field
from ..field.field_type import FieldType


class Frame:
    """One output state containing fields evaluated at one frame value."""

    def __init__(
        self,
        id: int = 0,
        *,
        value: float | None = None,
        name: str | None = None,
    ) -> None:
        self.id = int(id)
        self.value = None if value is None else float(value)
        self.name = name
        self.fields: dict[str, Field] = {}

    def add(self, field: Field) -> Field:
        """Add or replace a field by its native semantic name."""

        self.fields[field.name] = field
        return field

    def field(self, key: str | FieldType) -> Field:
        """Return a field by exact name or canonical semantic ``FieldType``."""

        if isinstance(key, str):
            try:
                return self.fields[key]
            except KeyError:
                semantic = FieldType.from_name(key)
                if semantic is not FieldType.UNKNOWN:
                    matches = [
                        field
                        for field in self.fields.values()
                        if field.type is semantic
                    ]
                    if len(matches) == 1:
                        return matches[0]
                raise

        matches = [
            field
            for field in self.fields.values()
            if field.type is key
        ]
        if len(matches) == 1:
            return matches[0]
        if not matches:
            raise KeyError(f"frame has no field of type {key.value}")
        raise KeyError(
            f"frame has multiple fields of type {key.value}; use the exact name"
        )

    def __getitem__(self, key: str | FieldType) -> Field:
        return self.field(key)

    def __len__(self) -> int:
        return len(self.fields)
