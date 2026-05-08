"""Model consistency checks for object ownership and references."""

from __future__ import annotations

from femaster_api.model.constraints import ConnectorConstraint, CouplingConstraint, RBMConstraint, TieConstraint
from femaster_api.model.elements import Element
from femaster_api.model.loads import InertialLoad, NodalForce, PressureLoad, SurfaceTraction, ThermalLoad, VolumeLoad
from femaster_api.model.nodes import Node
from femaster_api.model.sections import BeamSection, ShellSection, SolidSection, TrussSection
from femaster_api.model.sets import ElementSet, EntitySet, NodeSet, SurfaceSet
from femaster_api.model.surfaces import SurfaceDefinition
from femaster_api.validation.diagnostics import Diagnostics


def validate_model(model) -> Diagnostics:
    diagnostics = Diagnostics()

    if len(model.nodes) == 0:
        diagnostics.error("model has no nodes", location="nodes")
    if len(model.elements) == 0:
        diagnostics.warning("model has no elements", location="elements")

    for element in model.elements:
        for node in element.nodes:
            if node not in model.nodes:
                diagnostics.error("element references a node outside this model", location="elements")

    for item in model.node_sets:
        for member in item:
            if member not in model.nodes:
                diagnostics.error(f"node set {item.name!r} references a node outside this model", location="node_sets")

    for item in model.element_sets:
        for member in item:
            if member not in model.elements:
                diagnostics.error(f"element set {item.name!r} references an element outside this model", location="element_sets")

    for item in model.surface_sets:
        for member in item:
            if member not in model.surfaces:
                diagnostics.error(f"surface set {item.name!r} references a surface outside this model", location="surface_sets")

    for surface in model.surfaces:
        _check_target(model, diagnostics, surface.target, location="surfaces")

    for section in model.sections:
        if section.material and not model.materials.has(section.material):
            diagnostics.error(f"section {section.name!r} references material outside this model", location="sections")
        if section.element_set not in model.element_sets:
            diagnostics.error(f"section {section.name!r} references element set outside this model", location="sections")
        if isinstance(section, BeamSection) and not model.sections.has_profile(section.profile):
            diagnostics.error(f"beam section {section.name!r} references profile outside this model", location="sections")
        if isinstance(section, (SolidSection, ShellSection)) and section.orientation and not model.orientations.has(section.orientation):
            diagnostics.error(f"section {section.name!r} references orientation outside this model", location="sections")
        if isinstance(section, TrussSection) and section.area <= 0:
            diagnostics.error(f"truss section {section.name!r} area must be positive", location="sections")

    for load in model.loads:
        _check_load(model, diagnostics, load)

    for support in model.supports:
        if support.orientation and not model.orientations.has(support.orientation):
            diagnostics.error("support references orientation outside this model", location="supports")
        _check_target(model, diagnostics, support.target, location="supports")

    for collector in model.support_collectors:
        for support in collector.supports:
            _check_support(model, diagnostics, support, location="support_collectors")

    for collector in model.load_collectors:
        for load in collector.loads:
            _check_load(model, diagnostics, load)

    for loadcase in model.loadcases:
        for collector in getattr(loadcase, "loads", ()):
            if not model.load_collectors.has(collector):
                diagnostics.error(
                    f"loadcase {loadcase.name!r} references load collector outside this model",
                    location="loadcases",
                )
        for collector in getattr(loadcase, "supports", ()):
            if not model.support_collectors.has(collector):
                diagnostics.error(
                    f"loadcase {loadcase.name!r} references support collector outside this model",
                    location="loadcases",
                )
        density = getattr(loadcase, "density", None)
        orientation = getattr(loadcase, "orientation", None)
        initial_velocity = getattr(loadcase, "initial_velocity", None)
        if density and not model.fields.has(density):
            diagnostics.error(f"loadcase {loadcase.name!r} references density field outside this model", location="loadcases")
        if orientation and not model.fields.has(orientation):
            diagnostics.error(f"loadcase {loadcase.name!r} references orientation field outside this model", location="loadcases")
        if initial_velocity and not model.fields.has(initial_velocity):
            diagnostics.error(
                f"loadcase {loadcase.name!r} references initial velocity field outside this model",
                location="loadcases",
            )

    for field in model.fields:
        for key in field.values:
            _check_field_key(model, diagnostics, key, location="fields")

    for constraint in model.constraints:
        if isinstance(constraint, RBMConstraint) and constraint.element_set not in model.element_sets:
            diagnostics.error("RBM constraint references element set outside this model", location="constraints")
        if isinstance(constraint, ConnectorConstraint) and not model.orientations.has(constraint.coordinate_system):
            diagnostics.error("connector references coordinate system outside this model", location="constraints")
        if isinstance(constraint, CouplingConstraint):
            if constraint.master not in model.node_sets:
                diagnostics.error("coupling references master node set outside this model", location="constraints")
            if not _known_set(model, constraint.slave):
                diagnostics.error("coupling references slave set outside this model", location="constraints")
        if isinstance(constraint, TieConstraint):
            if constraint.master not in model.surface_sets:
                diagnostics.warning(f"tie master {constraint.master.name!r} is not a known surface set", location="constraints")
            if constraint.slave not in model.surface_sets:
                diagnostics.warning(f"tie slave {constraint.slave.name!r} is not a known surface set", location="constraints")

    return diagnostics


def _check_load(model, diagnostics: Diagnostics, load) -> None:
    _check_target(model, diagnostics, getattr(load, "target", None), location="loads")
    if isinstance(load, NodalForce):
        if load.orientation and not model.orientations.has(load.orientation):
            diagnostics.error("load references orientation outside this model", location="loads")
        if load.amplitude and not model.loads.has_amplitude(load.amplitude):
            diagnostics.error("load references amplitude outside this model", location="loads")
    elif isinstance(load, (SurfaceTraction, VolumeLoad)):
        if load.orientation and not model.orientations.has(load.orientation):
            diagnostics.error("load references orientation outside this model", location="loads")
        if load.amplitude and not model.loads.has_amplitude(load.amplitude):
            diagnostics.error("load references amplitude outside this model", location="loads")
    elif isinstance(load, PressureLoad):
        if load.amplitude and not model.loads.has_amplitude(load.amplitude):
            diagnostics.error("load references amplitude outside this model", location="loads")
    elif isinstance(load, ThermalLoad):
        if not model.fields.has(load.temperature_field):
            diagnostics.error("thermal load references field outside this model", location="loads")
    elif isinstance(load, InertialLoad):
        if load.target not in model.element_sets:
            diagnostics.error("inertial load references element set outside this model", location="loads")


def _check_support(model, diagnostics: Diagnostics, support, *, location: str) -> None:
    if support.orientation and not model.orientations.has(support.orientation):
        diagnostics.error("support references orientation outside this model", location=location)
    _check_target(model, diagnostics, support.target, location=location)


def _check_target(model, diagnostics: Diagnostics, target: object, *, location: str) -> None:
    if target is None:
        return
    if isinstance(target, Node) and target not in model.nodes:
        diagnostics.error("target references a node outside this model", location=location)
    elif isinstance(target, Element) and target not in model.elements:
        diagnostics.error("target references an element outside this model", location=location)
    elif isinstance(target, NodeSet) and target not in model.node_sets:
        diagnostics.error("target references a node set outside this model", location=location)
    elif isinstance(target, ElementSet) and target not in model.element_sets:
        diagnostics.error("target references an element set outside this model", location=location)
    elif isinstance(target, SurfaceDefinition) and target not in model.surfaces:
        diagnostics.error("target references a surface outside this model", location=location)
    elif isinstance(target, SurfaceSet) and target not in model.surface_sets:
        diagnostics.error("target references a surface set outside this model", location=location)


def _check_field_key(model, diagnostics: Diagnostics, key: object, *, location: str) -> None:
    if isinstance(key, tuple):
        for item in key:
            _check_field_key(model, diagnostics, item, location=location)
    else:
        _check_target(model, diagnostics, key, location=location)


def _known_set(model, entity_set: EntitySet) -> bool:
    if isinstance(entity_set, NodeSet):
        return entity_set in model.node_sets
    if isinstance(entity_set, ElementSet):
        return entity_set in model.element_sets
    if isinstance(entity_set, SurfaceSet):
        return entity_set in model.surface_sets
    return False
