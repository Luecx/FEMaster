"""Central field definitions shared by model input and result importers.

Field domains describe where rows live. Field types describe semantic meaning
independently of the result format that carried the values. This prevents RES,
FRD and future importers from maintaining separate name conventions.
"""

from __future__ import annotations

from enum import Enum
from typing import Iterable

from ._format import block, csv, keyword
from .repository import NamedObject, NamedRepository


class FieldDomain(Enum):
    """Physical storage domains implemented by FEMaster ModelData."""

    UNKNOWN       = "UNKNOWN"
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
    BUCKLING_MODE              = "BUCKLING_MODE"
    EIGENVALUE                 = "EIGENVALUE"
    PARTICIPATION              = "PARTICIPATION"
    VELOCITY                   = "VELOCITY"
    ACCELERATION               = "ACCELERATION"
    REACTION_FORCE             = "REACTION_FORCE"
    EXTERNAL_FORCE             = "EXTERNAL_FORCE"
    INTERNAL_FORCE             = "INTERNAL_FORCE"
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
        """Resolve current FEMaster, RES and FRD aliases to one canonical type.

        Transient, modal and buckling writers currently append frame/mode numbers
        to several field names. Those suffixes deliberately do not create new
        semantic field types: ``MODE_SHAPE_3`` still resolves to MODE_SHAPE.
        """

        normalized = "".join(character for character in name.upper() if character.isalnum())

        prefixes = (
            ("MODESHAPE", cls.MODE_SHAPE),
            ("BUCKLINGMODE", cls.BUCKLING_MODE),
            ("PARTICIPATION", cls.PARTICIPATION),
            ("DISPLACEMENT", cls.DISPLACEMENT),
            ("VELOCITY", cls.VELOCITY),
            ("ACCELERATION", cls.ACCELERATION),
        )
        for prefix, field_type in prefixes:
            if normalized.startswith(prefix):
                return field_type

        aliases = {
            "POSITION": cls.POSITION,
            "DISP": cls.DISPLACEMENT,
            "U": cls.DISPLACEMENT,
            "MODE": cls.MODE_SHAPE,
            "EIGENVALUE": cls.EIGENVALUE,
            "EIGENVALUES": cls.EIGENVALUE,
            "VELO": cls.VELOCITY,
            "ACCE": cls.ACCELERATION,
            "FORC": cls.REACTION_FORCE,
            "REACTIONFORCE": cls.REACTION_FORCE,
            "REACTIONFORCES": cls.REACTION_FORCE,
            "EXTFORC": cls.EXTERNAL_FORCE,
            "EXTERNALFORCE": cls.EXTERNAL_FORCE,
            "EXTERNALFORCES": cls.EXTERNAL_FORCE,
            "INTFORC": cls.INTERNAL_FORCE,
            "INTERNALFORCE": cls.INTERNAL_FORCE,
            "INTERNALFORCES": cls.INTERNAL_FORCE,
            "E": cls.STRAIN,
            "STRAIN": cls.STRAIN,
            "S": cls.STRESS,
            "STRESS": cls.STRESS,
            "STOP": cls.STRESS_TOP,
            "STRESSTOP": cls.STRESS_TOP,
            "SBOT": cls.STRESS_BOTTOM,
            "STRESSBOTTOM": cls.STRESS_BOTTOM,
            "SHR": cls.SHELL_RESULTANTS,
            "SHELLRESULTANTS": cls.SHELL_RESULTANTS,
            "SF": cls.SECTION_FORCE,
            "SECTIONFORCE": cls.SECTION_FORCE,
            "SHEAR": cls.SHEAR_FLOW,
            "SHEARFLOW": cls.SHEAR_FLOW,
            "PEEQ": cls.EQUIVALENT_PLASTIC_STRAIN,
        }
        if normalized in aliases:
            return aliases[normalized]

        for item in cls:
            if normalized in {
                "".join(character for character in item.name.upper() if character.isalnum()),
                "".join(character for character in item.value.upper() if character.isalnum()),
            }:
                return item
        return cls.UNKNOWN


FieldIndex = int | str
FieldKey = FieldIndex | tuple[FieldIndex, ...]


class Field(NamedObject):
    """Sparse field with explicit domain, semantics and component names.

    Values are keyed by their semantic model identifier. Node and element fields
    may use either bare local ids (``17``) or instance-qualified RES identifiers
    (``"bolt.17"``). Element-location fields use tuples such as
    ``("bolt.17", local_index)``.
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
        """Return the number of numerical components in one field row."""

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

    def export(self) -> str:
        """Export this field as a FIELD block."""

        lines = [
            keyword(
                "FIELD",
                NAME=self.name,
                TYPE=self.domain.value,
                COLS=self.cols,
                FILL="ZERO",
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

    def export(self) -> str:
        return "\n\n".join(field.export() for field in self)
