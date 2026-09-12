"""Canonical semantic field types."""

from enum import Enum


class FieldType(Enum):
    """Format-independent semantic field types."""

    UNKNOWN = "UNKNOWN"
    POSITION = "POSITION"
    DISPLACEMENT = "DISPLACEMENT"
    MODE_SHAPE = "MODE_SHAPE"
    BUCKLING_MODE = "BUCKLING_MODE"
    EIGENVALUE = "EIGENVALUE"
    PARTICIPATION = "PARTICIPATION"
    VELOCITY = "VELOCITY"
    ACCELERATION = "ACCELERATION"
    REACTION_FORCE = "REACTION_FORCE"
    EXTERNAL_FORCE = "EXTERNAL_FORCE"
    INTERNAL_FORCE = "INTERNAL_FORCE"
    TEMPERATURE = "TEMPERATURE"
    HEAT_FLUX = "HEAT_FLUX"
    STRESS = "STRESS"
    STRAIN = "STRAIN"
    GREEN_LAGRANGE_STRAIN = "GREEN_LAGRANGE_STRAIN"
    LOGARITHMIC_STRAIN = "LOGARITHMIC_STRAIN"
    PLASTIC_STRAIN = "PLASTIC_STRAIN"
    EQUIVALENT_PLASTIC_STRAIN = "PEEQ"
    STRESS_TOP = "STRESS_TOP"
    STRESS_BOTTOM = "STRESS_BOTTOM"
    SHELL_RESULTANTS = "SHELL_RESULTANTS"
    SECTION_FORCE = "SECTION_FORCE"
    SHEAR_FLOW = "SHEAR_FLOW"
    COMPLIANCE = "COMPLIANCE"
    VOLUME = "VOLUME"
    DENSITY = "DENSITY"
    MATERIAL_ORIENTATION = "MATERIAL_ORIENTATION"

    @classmethod
    def from_name(cls, name: str) -> "FieldType":
        normalized = "".join(
            character for character in name.upper() if character.isalnum()
        )

        for prefix, field_type in (
            ("MODESHAPE", cls.MODE_SHAPE),
            ("BUCKLINGMODE", cls.BUCKLING_MODE),
            ("PARTICIPATION", cls.PARTICIPATION),
            ("DISPLACEMENT", cls.DISPLACEMENT),
            ("VELOCITY", cls.VELOCITY),
            ("ACCELERATION", cls.ACCELERATION),
        ):
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
