"""Geometry and set commands for FEMaster input decks."""

from __future__ import annotations

from dataclasses import replace

from femaster_api.model import Element, ElementSet, ElementTopology, EntitySet, EntityType, Node, NodeSet, SurfaceSet
from femaster_api.model.surfaces import SurfaceDefinition

from .utils import FEMasterInputError, Header, parse_csv, side_token, truthy


ELEMENT_TYPES: dict[str, ElementTopology] = {
    "C3D4": ElementTopology.TET4,
    "C3D5": ElementTopology.PYRAMID5,
    "C3D6": ElementTopology.WEDGE6,
    "C3D8": ElementTopology.HEX8,
    "C3D10": ElementTopology.TET10,
    "C3D15": ElementTopology.WEDGE15,
    "C3D20": ElementTopology.HEX20,
    "C3D20R": ElementTopology.HEX20R,
    "T3": ElementTopology.TRUSS2,
    "B33": ElementTopology.BEAM2,
    "S3": ElementTopology.TRI3,
    "C2D3": ElementTopology.TRI3,
    "S4": ElementTopology.QUAD4,
    "C2D4": ElementTopology.QUAD4,
    "MITC4": ElementTopology.MITC4,
    "MITC4FRT": ElementTopology.MITC4,
    "QSPT": ElementTopology.QSPT,
    "S6": ElementTopology.TRI6,
    "C2D6": ElementTopology.TRI6,
    "S8": ElementTopology.QUAD8,
    "C2D8": ElementTopology.QUAD8,
}


def cmd_heading(parser, header: Header) -> None:
    for _ in parser.consume_data_lines():
        pass


def cmd_model(parser, header: Header) -> None:
    if "NAME" in header.params:
        parser.model.name = header.params["NAME"]


def cmd_node(parser, header: Header) -> None:
    node_set_name = header.params.get("NSET")
    for line in parser.consume_data_lines():
        tokens = parse_csv(line)
        if len(tokens) < 4:
            raise FEMasterInputError("*NODE rows require id,x,y,z")
        deck_id = int(tokens[0])
        if deck_id in parser.node_by_deck_id:
            raise FEMasterInputError(f"duplicate node id {deck_id}")
        node = parser.model.nodes.add(Node(float(tokens[1]), float(tokens[2]), float(tokens[3])))
        parser.node_by_deck_id[deck_id] = node
        _add_to_set(parser, EntityType.NODE, "NALL", node)
        if node_set_name:
            _add_to_set(parser, EntityType.NODE, node_set_name, node)


def cmd_element(parser, header: Header) -> None:
    raw_type = header.params.get("TYPE")
    if not raw_type:
        raise FEMasterInputError("*ELEMENT requires TYPE")
    try:
        topology = ELEMENT_TYPES[raw_type.upper()]
    except KeyError as exc:
        raise FEMasterInputError(f"unsupported element type {raw_type!r}") from exc
    element_set_name = header.params.get("ELSET")

    for line in parser.consume_data_lines():
        tokens = [token for token in parse_csv(line) if token]
        if len(tokens) < 2:
            continue
        deck_id = int(tokens[0])
        if deck_id in parser.element_by_deck_id:
            raise FEMasterInputError(f"duplicate element id {deck_id}")
        try:
            nodes = tuple(parser.node_by_deck_id[int(token)] for token in tokens[1:])
        except KeyError as exc:
            raise FEMasterInputError(f"element {deck_id} references unknown node {exc.args[0]}") from exc
        element = parser.model.elements.add(Element(topology, nodes))
        parser.element_by_deck_id[deck_id] = element
        _add_to_set(parser, EntityType.ELEMENT, "EALL", element)
        if element_set_name:
            _add_to_set(parser, EntityType.ELEMENT, element_set_name, element)


def cmd_nset(parser, header: Header) -> None:
    name = header.params.get("NSET") or header.params.get("NAME")
    if not name:
        raise FEMasterInputError("*NSET requires NSET or NAME")
    for deck_id in _ids_from_lines(parser, truthy(header.params.get("GENERATE"))):
        try:
            node = parser.node_by_deck_id[deck_id]
        except KeyError as exc:
            raise FEMasterInputError(f"*NSET {name} references unknown node {deck_id}") from exc
        _add_to_set(parser, EntityType.NODE, name, node, generated=truthy(header.params.get("GENERATE")))


def cmd_elset(parser, header: Header) -> None:
    name = header.params.get("ELSET") or header.params.get("NAME")
    if not name:
        raise FEMasterInputError("*ELSET requires ELSET or NAME")
    for deck_id in _ids_from_lines(parser, truthy(header.params.get("GENERATE"))):
        try:
            element = parser.element_by_deck_id[deck_id]
        except KeyError as exc:
            raise FEMasterInputError(f"*ELSET {name} references unknown element {deck_id}") from exc
        _add_to_set(parser, EntityType.ELEMENT, name, element, generated=truthy(header.params.get("GENERATE")))


def cmd_surface(parser, header: Header) -> None:
    name = header.params.get("NAME") or header.params.get("SFSET")
    if not name:
        raise FEMasterInputError("*SURFACE requires NAME or SFSET")
    for line in parser.consume_data_lines():
        tokens = [token for token in parse_csv(line) if token]
        if len(tokens) == 3:
            sid = int(tokens[0])
            target = _surface_target(parser, tokens[1])
            side = side_token(tokens[2])
        elif len(tokens) == 2:
            sid = None
            target = _surface_target(parser, tokens[0])
            side = side_token(tokens[1])
        else:
            raise FEMasterInputError("*SURFACE rows require target,side or id,target,side")
        surface = parser.model.surfaces.add(SurfaceDefinition(name, target, side, sid))
        if sid is not None:
            parser.surface_by_deck_id[sid] = surface
        _add_to_set(parser, EntityType.SURFACE, name, surface)


def cmd_sfset(parser, header: Header) -> None:
    name = header.params.get("SFSET") or header.params.get("NAME")
    if not name:
        raise FEMasterInputError("*SFSET requires SFSET or NAME")
    for deck_id in _ids_from_lines(parser, truthy(header.params.get("GENERATE"))):
        try:
            surface = parser.surface_by_deck_id[deck_id]
        except KeyError as exc:
            raise FEMasterInputError(f"*SFSET {name} references unknown surface {deck_id}") from exc
        _add_to_set(parser, EntityType.SURFACE, name, surface, generated=truthy(header.params.get("GENERATE")))


def _surface_target(parser, token: str):
    try:
        return parser.model.sets.get(EntityType.ELEMENT, token)
    except Exception:
        pass
    try:
        return parser.element_by_deck_id[int(token)]
    except (KeyError, ValueError) as exc:
        raise FEMasterInputError(f"surface target {token!r} is neither an ELSET nor an element id") from exc


def _ids_from_lines(parser, generated: bool):
    for line in parser.consume_data_lines():
        tokens = [token for token in parse_csv(line) if token]
        if generated:
            if len(tokens) < 2:
                raise FEMasterInputError("GENERATE set rows require start,end[,step]")
            start = int(tokens[0])
            end = int(tokens[1])
            step = int(tokens[2]) if len(tokens) > 2 else 1
            if step == 0:
                raise FEMasterInputError("GENERATE step must not be zero")
            current = start
            while (step > 0 and current <= end) or (step < 0 and current >= end):
                yield current
                current += step
        else:
            for token in tokens:
                yield int(token)


def _set_class(entity_type: EntityType):
    return {
        EntityType.NODE: NodeSet,
        EntityType.ELEMENT: ElementSet,
        EntityType.SURFACE: SurfaceSet,
    }[entity_type]


def _add_to_set(parser, entity_type: EntityType, name: str, member: object, *, generated: bool = False) -> EntitySet:
    try:
        current = parser.model.sets.get(entity_type, name)
    except KeyError:
        current = parser.model.sets.add(_set_class(entity_type)(name, ()))
    updated = replace(current, members=(*current.members, member), generated=current.generated or generated)
    return parser.model.sets.add(updated)


GEOMETRY_COMMANDS = {
    "HEADING": cmd_heading,
    "SUMMARY": cmd_heading,
    "MODEL": cmd_model,
    "ASSEMBLY": lambda parser, header: None,
    "ENDASSEMBLY": lambda parser, header: None,
    "PART": lambda parser, header: None,
    "ENDPART": lambda parser, header: None,
    "NODE": cmd_node,
    "NODECOUNT": lambda parser, header: [line for line in parser.consume_data_lines()],
    "ELEMENT": cmd_element,
    "ELEMENTCOUNT": lambda parser, header: [line for line in parser.consume_data_lines()],
    "NSET": cmd_nset,
    "ELSET": cmd_elset,
    "SURFACE": cmd_surface,
    "SURFACECOUNT": lambda parser, header: [line for line in parser.consume_data_lines()],
    "SFSET": cmd_sfset,
}
