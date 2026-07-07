"""Load, support, amplitude, and constraint commands."""

from __future__ import annotations

from femaster_api.model import (
    Amplitude,
    AmplitudeInterpolation,
    ConnectorConstraint,
    ConnectorType,
    CouplingConstraint,
    CouplingType,
    InertialLoad,
    LoadCollector,
    NodalForce,
    PressureLoad,
    RBMConstraint,
    Support,
    SupportCollector,
    SurfaceTraction,
    ThermalLoad,
    TieConstraint,
    VolumeLoad,
)
from femaster_api.model.sets import EntityType

from .utils import FEMasterInputError, Header, falsey, numbers, parse_csv, truthy


def cmd_amplitude(parser, header: Header) -> None:
    name = header.params.get("NAME")
    if not name:
        raise FEMasterInputError("*AMPLITUDE requires NAME")
    interpolation = AmplitudeInterpolation[header.params.get("TYPE", "LINEAR").upper()]
    values = [value for line in parser.consume_data_lines() for value in numbers(line)]
    if len(values) % 2 != 0:
        raise FEMasterInputError("*AMPLITUDE requires time,value pairs")
    points = tuple((values[index], values[index + 1]) for index in range(0, len(values), 2))
    parser.model.loads.add_amplitude(Amplitude(name, points, interpolation))


def cmd_support(parser, header: Header) -> None:
    name = header.params.get("SUPPORT_COLLECTOR")
    if not name:
        raise FEMasterInputError("*SUPPORT requires SUPPORT_COLLECTOR")
    orientation = _orientation(parser, header)
    collector = _support_collector(parser, name)
    for line in parser.consume_data_lines():
        tokens = parse_csv(line)
        if not tokens:
            continue
        target = _node_target(parser, tokens[0])
        values = tuple(_support_value(token) for token in (tokens[1:7] + [""] * 6)[:6])
        collector = collector.add(Support(target, values, orientation))
    parser.model.support_collectors.add(collector)


def cmd_cload(parser, header: Header) -> None:
    collector = _load_collector(parser, header)
    orientation = _orientation(parser, header)
    amplitude = _amplitude(parser, header)
    for line in parser.consume_data_lines():
        tokens = parse_csv(line)
        if not tokens:
            continue
        target = _node_target(parser, tokens[0])
        collector = collector.add(NodalForce(target, _vector(tokens[1:], 6), orientation, amplitude))
    parser.model.load_collectors.add(collector)


def cmd_dload(parser, header: Header) -> None:
    collector = _load_collector(parser, header)
    orientation = _orientation(parser, header)
    amplitude = _amplitude(parser, header)
    for line in parser.consume_data_lines():
        tokens = parse_csv(line)
        if len(tokens) < 2:
            continue
        collector = collector.add(SurfaceTraction(_surface_target(parser, tokens[0]), _vector(tokens[1:], 3), orientation, amplitude))
    parser.model.load_collectors.add(collector)


def cmd_pload(parser, header: Header) -> None:
    collector = _load_collector(parser, header)
    amplitude = _amplitude(parser, header)
    for line in parser.consume_data_lines():
        tokens = parse_csv(line)
        if len(tokens) < 2:
            continue
        collector = collector.add(PressureLoad(_surface_target(parser, tokens[0]), float(tokens[1]), amplitude))
    parser.model.load_collectors.add(collector)


def cmd_vload(parser, header: Header) -> None:
    collector = _load_collector(parser, header)
    orientation = _orientation(parser, header)
    amplitude = _amplitude(parser, header)
    for line in parser.consume_data_lines():
        tokens = parse_csv(line)
        if len(tokens) < 2:
            continue
        collector = collector.add(VolumeLoad(_element_target(parser, tokens[0]), _vector(tokens[1:], 3), orientation, amplitude))
    parser.model.load_collectors.add(collector)


def cmd_tload(parser, header: Header) -> None:
    collector = _load_collector(parser, header)
    field_name = header.params.get("TEMPERATUREFIELD") or header.params.get("FIELD")
    if not field_name:
        raise FEMasterInputError("*TLOAD requires TEMPERATUREFIELD")
    ref = header.params.get("REFERENCETEMPERATURE")
    if ref is None:
        raise FEMasterInputError("*TLOAD requires REFERENCETEMPERATURE")
    _reject_data(parser, "*TLOAD")
    collector = collector.add(ThermalLoad(parser.model.fields.get(field_name), float(ref)))
    parser.model.load_collectors.add(collector)


def cmd_inertial_load(parser, header: Header) -> None:
    collector = _load_collector(parser, header)
    consider = truthy(header.params.get("CONSIDER_POINT_MASSES"))
    for line in parser.consume_data_lines():
        tokens = parse_csv(line)
        if len(tokens) < 13:
            raise FEMasterInputError("*INERTIALOAD rows require target plus 12 numeric values")
        target = parser.model.sets.get(EntityType.ELEMENT, tokens[0])
        values = [float(token) for token in tokens[1:13]]
        collector = collector.add(InertialLoad(target, tuple(values[0:3]), tuple(values[3:6]), tuple(values[6:9]), tuple(values[9:12]), consider))
    parser.model.load_collectors.add(collector)


def cmd_rbm(parser, header: Header) -> None:
    name = header.params.get("ELSET") or "EALL"
    _reject_data(parser, "*RBM")
    parser.model.constraints.add(RBMConstraint(parser.model.sets.get(EntityType.ELEMENT, name)))


def cmd_coupling(parser, header: Header) -> None:
    master = header.params.get("MASTER")
    slave = header.params.get("SLAVE")
    sfset = header.params.get("SFSET")
    if not master or bool(slave) == bool(sfset):
        raise FEMasterInputError("*COUPLING requires MASTER and exactly one of SLAVE/SFSET")
    dof_lines = list(parser.consume_data_lines())
    if not dof_lines:
        raise FEMasterInputError("*COUPLING requires a DOF line")
    dofs = tuple(bool(int(float(token))) for token in parse_csv(dof_lines[0])[:6])
    if len(dofs) != 6:
        raise FEMasterInputError("*COUPLING DOF line requires six values")
    coupling_type = CouplingType[header.params.get("TYPE", "KINEMATIC").upper().replace("DISTRIBUTING", "STRUCTURAL")]
    master_set = parser.model.sets.get(EntityType.NODE, master)
    slave_set = parser.model.sets.get(EntityType.SURFACE if sfset else EntityType.NODE, sfset or slave)
    parser.model.constraints.add(CouplingConstraint(master_set, slave_set, dofs, coupling_type))


def cmd_connector(parser, header: Header) -> None:
    ctype = header.params.get("TYPE")
    nset1 = header.params.get("NSET1")
    nset2 = header.params.get("NSET2")
    coords = header.params.get("COORDINATESYSTEM")
    if not (ctype and nset1 and nset2 and coords):
        raise FEMasterInputError("*CONNECTOR requires TYPE,NSET1,NSET2,COORDINATESYSTEM")
    parser.model.constraints.add(
        ConnectorConstraint(
            ConnectorType[ctype.upper()],
            parser.model.sets.get(EntityType.NODE, nset1),
            parser.model.sets.get(EntityType.NODE, nset2),
            parser.model.orientations.get(coords),
        )
    )


def cmd_tie(parser, header: Header) -> None:
    master = header.params.get("MASTER")
    slave = header.params.get("SLAVE")
    if not (master and slave):
        raise FEMasterInputError("*TIE requires MASTER and SLAVE")
    _reject_data(parser, "*TIE")
    parser.model.constraints.add(
        TieConstraint(
            parser.model.sets.get(EntityType.SURFACE, master),
            _tie_slave(parser, slave),
            float(header.params.get("DISTANCE", "0.0")),
            truthy(header.params.get("ADJUST")) and not falsey(header.params.get("ADJUST")),
        )
    )


def _load_collector(parser, header: Header) -> LoadCollector:
    name = header.params.get("LOAD_COLLECTOR")
    if not name:
        raise FEMasterInputError(f"*{header.keyword} requires LOAD_COLLECTOR")
    try:
        return parser.model.load_collectors.get(name)
    except KeyError:
        return parser.model.load_collectors.add(LoadCollector(name))


def _support_collector(parser, name: str) -> SupportCollector:
    try:
        return parser.model.support_collectors.get(name)
    except KeyError:
        return parser.model.support_collectors.add(SupportCollector(name))


def _orientation(parser, header: Header):
    name = header.params.get("ORIENTATION")
    return None if not name else parser.model.orientations.get(name)


def _amplitude(parser, header: Header):
    name = header.params.get("AMPLITUDE")
    if not name:
        return None
    for amplitude in parser.model.loads.amplitudes():
        if amplitude.name == name:
            return amplitude
    raise FEMasterInputError(f"unknown amplitude {name!r}")


def _node_target(parser, token: str):
    try:
        return parser.model.sets.get(EntityType.NODE, token)
    except Exception:
        pass
    try:
        return parser.node_by_deck_id[int(token)]
    except (KeyError, ValueError) as exc:
        raise FEMasterInputError(f"target {token!r} is neither an NSET nor a node id") from exc


def _element_target(parser, token: str):
    try:
        return parser.model.sets.get(EntityType.ELEMENT, token)
    except Exception:
        pass
    try:
        return parser.element_by_deck_id[int(token)]
    except (KeyError, ValueError) as exc:
        raise FEMasterInputError(f"target {token!r} is neither an ELSET nor an element id") from exc


def _surface_target(parser, token: str):
    try:
        return parser.model.sets.get(EntityType.SURFACE, token)
    except Exception:
        pass
    try:
        surface = parser.surface_by_deck_id[int(token)]
    except (KeyError, ValueError) as exc:
        raise FEMasterInputError(f"target {token!r} is neither an SFSET nor a surface id") from exc
    return surface


def _tie_slave(parser, token: str):
    try:
        return parser.model.sets.get(EntityType.NODE, token)
    except Exception:
        return parser.model.sets.get(EntityType.SURFACE, token)


def _vector(tokens: list[str], width: int) -> tuple[float, ...]:
    values = [0.0 if token == "" else float(token) for token in tokens[:width]]
    values.extend([0.0] * (width - len(values)))
    return tuple(values)


def _support_value(token: str) -> float | None:
    return None if token == "" else float(token)


def _reject_data(parser, command: str) -> None:
    for line in parser.consume_data_lines():
        if line.strip():
            raise FEMasterInputError(f"{command} does not take data lines")


LOAD_COMMANDS = {
    "AMPLITUDE": cmd_amplitude,
    "SUPPORT": cmd_support,
    "CLOAD": cmd_cload,
    "DLOAD": cmd_dload,
    "PLOAD": cmd_pload,
    "VLOAD": cmd_vload,
    "TLOAD": cmd_tload,
    "INERTIALOAD": cmd_inertial_load,
    "INERTIALLOAD": cmd_inertial_load,
    "RBM": cmd_rbm,
    "COUPLING": cmd_coupling,
    "CONNECTOR": cmd_connector,
    "TIE": cmd_tie,
}
