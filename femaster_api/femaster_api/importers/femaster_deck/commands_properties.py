"""Material, field, orientation, feature, and section commands."""

from __future__ import annotations

from femaster_api.model import (
    BeamSection,
    CylindricalOrientation,
    Field,
    FieldDomain,
    Material,
    PointMass,
    Profile,
    RectangularOrientation,
    ShellSection,
    SolidSection,
    TrussSection,
)
from femaster_api.model.materials import ABDElasticity, GeneralizedIsotropicElasticity, IsotropicElasticity, OrthotropicElasticity
from femaster_api.model.sets import EntityType

from .utils import FEMasterInputError, Header, numbers, parse_csv


def cmd_material(parser, header: Header) -> None:
    name = header.params.get("NAME") or header.params.get("MATERIAL")
    if not name:
        raise FEMasterInputError("*MATERIAL requires NAME")
    material = Material(name)
    parser.current_material = parser.model.materials.add(material)


def cmd_elastic(parser, header: Header) -> None:
    material = _active_material(parser)
    values = [value for line in parser.consume_data_lines() for value in numbers(line)]
    if not values:
        raise FEMasterInputError("*ELASTIC requires data")

    kind = header.params.get("TYPE", "ISOTROPIC").replace("_", "").upper()
    if kind in {"ISO", "ISOTROPIC"}:
        if len(values) < 2:
            raise FEMasterInputError("*ELASTIC, TYPE=ISOTROPIC requires E,nu")
        elasticity = IsotropicElasticity(values[0], values[1])
    elif kind in {"GENISO", "GENERALIZEDISOTROPIC"}:
        if len(values) < 3:
            raise FEMasterInputError("*ELASTIC, TYPE=GENISO requires E,nu,G")
        elasticity = GeneralizedIsotropicElasticity(values[0], values[1], values[2])
    elif kind in {"ORTHO", "ORTHOTROPIC"}:
        if len(values) < 9:
            raise FEMasterInputError("*ELASTIC, TYPE=ORTHOTROPIC requires 9 values")
        elasticity = OrthotropicElasticity(*values[:9])
    elif kind in {"ABD", "ABDELASTICITY"}:
        elasticity = ABDElasticity(tuple(values))
    else:
        raise FEMasterInputError(f"unsupported *ELASTIC TYPE={header.params.get('TYPE')!r}")

    parser.current_material = parser.model.materials.add(material.set_elasticity(elasticity))


def cmd_density(parser, header: Header) -> None:
    material = _active_material(parser)
    data = [value for line in parser.consume_data_lines() for value in numbers(line)]
    if not data:
        raise FEMasterInputError("*DENSITY requires one value")
    parser.current_material = parser.model.materials.add(material.set_density(data[0]))


def cmd_thermal_expansion(parser, header: Header) -> None:
    material = _active_material(parser)
    data = [value for line in parser.consume_data_lines() for value in numbers(line)]
    if not data:
        raise FEMasterInputError("*THERMALEXPANSION requires one value")
    parser.current_material = parser.model.materials.add(material.set_thermal_expansion(data[0]))


def cmd_orientation(parser, header: Header) -> None:
    name = header.params.get("NAME")
    if not name:
        raise FEMasterInputError("*ORIENTATION requires NAME")
    data = [value for line in parser.consume_data_lines() for value in numbers(line)]
    kind = header.params.get("TYPE", "RECTANGULAR").upper()
    if kind == "RECTANGULAR":
        if len(data) not in {3, 6, 9}:
            raise FEMasterInputError("rectangular *ORIENTATION requires 3, 6, or 9 values")
        orientation = RectangularOrientation(name, tuple(data[:3]), tuple(data[3:6]) if len(data) >= 6 else None, tuple(data[6:9]) if len(data) >= 9 else None)
    elif kind == "CYLINDRICAL":
        if len(data) < 9:
            raise FEMasterInputError("cylindrical *ORIENTATION requires 9 values")
        orientation = CylindricalOrientation(name, tuple(data[:3]), tuple(data[3:6]), tuple(data[6:9]))
    else:
        raise FEMasterInputError(f"unsupported *ORIENTATION TYPE={kind!r}")
    parser.model.orientations.add(orientation)


def cmd_profile(parser, header: Header) -> None:
    name = header.params.get("NAME") or header.params.get("PROFILE")
    if not name:
        raise FEMasterInputError("*PROFILE requires NAME")
    data = [value for line in parser.consume_data_lines() for value in numbers(line)]
    if len(data) < 4:
        raise FEMasterInputError("*PROFILE requires at least area,iy,iz,j")
    padded = data + [0.0] * (9 - len(data))
    parser.model.profiles.add(Profile(name, *padded[:9]))


def cmd_solid_section(parser, header: Header) -> None:
    name = _section_name(header, "SOLID")
    parser.model.sections.add(
        SolidSection(
            name,
            _material(parser, header),
            _element_set(parser, header),
            _orientation(parser, header),
        )
    )
    _reject_data(parser, "*SOLIDSECTION")


def cmd_shell_section(parser, header: Header) -> None:
    data = [value for line in parser.consume_data_lines() for value in numbers(line)]
    thickness = data[0] if data else float(header.params.get("THICKNESS", 1.0))
    parser.model.sections.add(
        ShellSection(
            _section_name(header, "SHELL"),
            _material(parser, header),
            _element_set(parser, header),
            thickness,
            _orientation(parser, header),
        )
    )


def cmd_beam_section(parser, header: Header) -> None:
    data = [value for line in parser.consume_data_lines() for value in numbers(line)]
    profile_name = header.params.get("PROFILE")
    if not profile_name:
        raise FEMasterInputError("*BEAMSECTION requires PROFILE")
    orientation = tuple(data[:3]) if len(data) >= 3 else None
    parser.model.sections.add(
        BeamSection(
            _section_name(header, "BEAM"),
            _material(parser, header),
            _element_set(parser, header),
            parser.model.profiles.get(profile_name),
            orientation,
        )
    )


def cmd_truss_section(parser, header: Header) -> None:
    data = [value for line in parser.consume_data_lines() for value in numbers(line)]
    area = data[0] if data else float(header.params.get("AREA", 0.0))
    parser.model.sections.add(TrussSection(_section_name(header, "TRUSS"), _material(parser, header), _element_set(parser, header), area))


def cmd_point_mass(parser, header: Header) -> None:
    nset = header.params.get("NSET")
    if not nset:
        raise FEMasterInputError("*POINTMASS requires NSET")
    data = [value for line in parser.consume_data_lines() for value in numbers(line)]
    data = data + [0.0] * (10 - len(data))
    parser.model.features.add(
        PointMass(
            parser.model.sets.get(EntityType.NODE, nset),
            data[0],
            tuple(data[1:4]),
            tuple(data[4:7]),
            tuple(data[7:10]),
        )
    )


def cmd_field(parser, header: Header) -> None:
    name = header.params.get("NAME")
    if not name:
        raise FEMasterInputError("*FIELD requires NAME")
    domain = _field_domain(header.params.get("TYPE") or header.params.get("DOMAIN"))
    cols = int(header.params.get("COLS", "0") or "0")
    rows: list[tuple[int | tuple[int, ...], tuple[float | None, ...]]] = []
    for line in parser.consume_data_lines():
        tokens = [token for token in parse_csv(line)]
        if not tokens:
            continue
        index_cols = _field_index_cols(domain)
        if len(tokens) <= index_cols:
            raise FEMasterInputError(f"*FIELD {name} row has no values")
        key = _field_key(parser, domain, tokens[:index_cols])
        key_value: int | tuple[int, ...] = key[0] if len(key) == 1 else key
        values = tuple(None if token == "" else float(token) for token in tokens[index_cols:])
        cols = cols or len(values)
        rows.append((key_value, values))
    if cols <= 0:
        raise FEMasterInputError(f"*FIELD {name} has no value columns")
    field = Field(name, domain, cols, fill=header.params.get("FILL", "ZERO"))
    for key, values in rows:
        field = field.set(key, values)
    parser.model.fields.add(field)


def cmd_ignored_property(parser, header: Header) -> None:
    for _ in parser.consume_data_lines():
        pass


def _active_material(parser) -> Material:
    if parser.current_material is None:
        raise FEMasterInputError("material property command without active *MATERIAL")
    return parser.current_material


def _section_name(header: Header, prefix: str) -> str:
    return header.params.get("NAME") or f"{prefix}_{header.params.get('ELSET', 'SECTION')}"


def _material(parser, header: Header) -> Material:
    name = header.params.get("MATERIAL") or header.params.get("MAT")
    if not name:
        raise FEMasterInputError(f"*{header.keyword} requires MATERIAL")
    return parser.model.materials.get(name)


def _element_set(parser, header: Header):
    name = header.params.get("ELSET")
    if not name:
        raise FEMasterInputError(f"*{header.keyword} requires ELSET")
    return parser.model.sets.get(EntityType.ELEMENT, name)


def _orientation(parser, header: Header):
    name = header.params.get("ORIENTATION")
    return None if not name else parser.model.orientations.get(name)


def _reject_data(parser, command: str) -> None:
    for line in parser.consume_data_lines():
        if line.strip():
            raise FEMasterInputError(f"{command} does not take data lines")


def _field_domain(raw: str | None) -> FieldDomain:
    if not raw:
        return FieldDomain.UNKNOWN
    normalized = raw.replace("_", "").upper()
    for domain in FieldDomain:
        if domain.name.replace("_", "").upper() == normalized or domain.value.replace("_", "").upper() == normalized:
            return domain
    if normalized == "IP":
        return FieldDomain.ELEMENT_IP
    raise FEMasterInputError(f"unsupported FIELD TYPE={raw!r}")


def _field_index_cols(domain: FieldDomain) -> int:
    if domain in {FieldDomain.NODE, FieldDomain.ELEMENT, FieldDomain.UNKNOWN}:
        return 1
    if domain in {FieldDomain.ELEMENT_NODAL, FieldDomain.ELEMENT_IP}:
        return 2
    return 1


def _field_key(parser, domain: FieldDomain, tokens: list[str]) -> tuple[int, ...]:
    if domain is FieldDomain.ELEMENT_NODAL:
        return (
            parser.element_by_deck_id[int(tokens[0])].id,
            parser.node_by_deck_id[int(tokens[1])].id,
        )
    if domain is FieldDomain.ELEMENT_IP:
        return (parser.element_by_deck_id[int(tokens[0])].id, int(tokens[1]))
    return (_model_field_index(parser, domain, int(tokens[0])),)


def _model_field_index(parser, domain: FieldDomain, deck_id: int) -> int:
    if domain is FieldDomain.NODE:
        return parser.node_by_deck_id[deck_id].id
    if domain in {FieldDomain.ELEMENT, FieldDomain.ELEMENT_NODAL, FieldDomain.ELEMENT_IP}:
        return parser.element_by_deck_id[deck_id].id
    return deck_id


PROPERTY_COMMANDS = {
    "MATERIAL": cmd_material,
    "ELASTIC": cmd_elastic,
    "DENSITY": cmd_density,
    "THERMALEXPANSION": cmd_thermal_expansion,
    "CONDUCTIVITY": cmd_ignored_property,
    "SPECIFICHEAT": cmd_ignored_property,
    "ORIENTATION": cmd_orientation,
    "PROFILE": cmd_profile,
    "SOLIDSECTION": cmd_solid_section,
    "SHELLSECTION": cmd_shell_section,
    "BEAMSECTION": cmd_beam_section,
    "TRUSSSECTION": cmd_truss_section,
    "POINTMASS": cmd_point_mass,
    "FIELD": cmd_field,
}
