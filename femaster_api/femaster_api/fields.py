"""Central field definitions shared by model input and result readers.

Field domains describe where rows live. Field types describe semantic meaning
independently of the result format that carried the values. This prevents RES,
FRD and future readers from maintaining separate name conventions.
"""

from __future__ import annotations

from enum import Enum
from typing import Iterable

from ._format import block, csv, keyword
from .repository import NamedObject, NamedRepository


class FieldDomain(Enum):
    """Physical storage domain of one FEMaster field."""

    UNKNOWN       = "UNKNOWN"
    GLOBAL        = "GLOBAL"
    NODE          = "NODE"
    ELEMENT       = "ELEMENT"
    ELEMENT_NODAL = "ELEMENT_NODAL"
    ELEMENT_IP    = "ELEMENT_IP"
    ELEMENT_MP    = "ELEMENT_MP"


class FieldType(Enum):
    """Canonical semantic field types used throughout the Python API."""

    UNKNOWN                    = "UNKNOWN"
    POSITION                   = "POSITION"
    DISPLACEMENT               = "DISPLACEMENT"
    MODE_SHAPE                 = "MODE_SHAPE"
    VELOCITY                   = "VELOCITY"
    ACCELERATION               = "ACCELERATION"
    REACTION_FORCE             = "REACTION_FORCE"
    TEMPERATURE                = "TEMPERATURE"
    HEAT_FLUX                  = "HEAT_FLUX"
    STRESS                     = "STRESS"
    STRAIN                     = "STRAIN"
    GREEN_LAGRANGE_STRAIN      = "GREEN_LAGRANGE_STRAIN"
    LOGARITHMIC_STRAIN         = "LOGARITHMIC_STRAIN"
    PLASTIC_STRAIN             = "PLASTIC_STRAIN"
    EQUIVALENT_PLASTIC_STRAIN = "PEEQ"
    STRESS_TOP                 = "STRESS_TOP"
    STRESS_BOTTOM              = "STRESS_BOTTOM"
    SHELL_RESULTANTS           = "SHELL_RESULTANTS"
    SECTION_FORCE              = "SECTION_FORCE"
    SHEAR_FLOW                 = "SHEAR_FLOW"
    COMPLIANCE                 = "COMPLIANCE"
    VOLUME                     = "VOLUME"
    DENSITY                    = "DENSITY"
    MATERIAL_ORIENTATION       = "MATERIAL_ORIENTATION"

    @classmethod
    def from_name(cls, name: str) -> "FieldType":
        """Resolve common FEMaster/FRD field names to one canonical type."""

        normalized = name.strip().upper().replace("-", "_").replace(" ", "_")
        aliases = {
            "DISP": cls.DISPLACEMENT,
            "U": cls.DISPLACEMENT,
            "MODE": cls.MODE_SHAPE,
            "MODE_SHAPE": cls.MODE_SHAPE,
            "E": cls.STRAIN,
            "S": cls.STRESS,
            "STRESS": cls.STRESS,
            "S_TOP": cls.STRESS_TOP,
            "S_BOT": cls.STRESS_BOTTOM,
            "SHR": cls.SHELL_RESULTANTS,
            "SF": cls.SECTION_FORCE,
            "SHEAR": cls.SHEAR_FLOW,
            "PEEQ": cls.EQUIVALENT_PLASTIC_STRAIN,
        }
        if normalized in aliases:
            return aliases[normalized]

        for item in cls:
            if normalized in {item.name, item.value}:
                return item
        return cls.UNKNOWN


FieldKey = int | tuple[int, ...] | str


class Field(NamedObject):
    """Sparse field with explicit domain, semantics and component names.

    Values are keyed by their semantic model identifier. Node and element fields
    normally use an integer key; element-nodal and integration-point fields may
    use tuples such as ``(element_id, local_index)``.
    """

    def __init__(
        self,
        name: str,
        domain: FieldDomain,
        components: Iterable[str] = (),
        *,
        type: FieldType | None = None,
    ) -> None:
        super().__init__(name)
        self.domain     = domain
        self.type       = type or FieldType.from_name(name)
        self.components = tuple(str(component) for component in components)
        self.values: dict[FieldKey, tuple[float, ...]] = {}

    @property
    def cols(self) -> int:
        """Return the number of value columns represented by each row."""

        if self.components:
            return len(self.components)
        if self.values:
            return len(next(iter(self.values.values())))
        return 0

    def set(self, key: FieldKey, values: Iterable[float]) -> "Field":
        """Set one field row and return this field for fluent construction."""

        row = tuple(float(value) for value in values)
        if self.components and len(row) != len(self.components):
            raise ValueError(
                f"field {self.name!r} expects {len(self.components)} values, got {len(row)}"
            )
        if self.values and len(row) != self.cols:
            raise ValueError(f"field {self.name!r} has inconsistent row width")
        self.values[key] = row
        return self

    def get(self, key: FieldKey) -> tuple[float, ...]:
        """Return one field row by semantic key."""

        return self.values[key]

    def to_femaster(self) -> str:
        """Return this field as a FEMaster FIELD block."""

        lines = [
            keyword(
                "FIELD",
                NAME=self.name,
                TYPE=self.domain.value,
                COLS=self.cols,
                FILL="NONE",
            )
        ]

        def sort_key(key: FieldKey) -> tuple[str, ...]:
            values = key if isinstance(key, tuple) else (key,)
            return tuple(str(value) for value in values)

        for key in sorted(self.values, key=sort_key):
            row_key = key if isinstance(key, tuple) else (key,)
            lines.append(csv((*row_key, *self.values[key])))

        return block(lines)

    def __getitem__(self, key: FieldKey) -> tuple[float, ...]:
        return self.get(key)

    def __len__(self) -> int:
        return len(self.values)


class FieldRepository(NamedRepository[Field]):
    """Named repository of global model fields."""

    def to_femaster(self) -> str:
        return "\n\n".join(field.to_femaster() for field in self)
