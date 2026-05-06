"""Model consistency checks."""

from __future__ import annotations

from femaster_api.model.constraints import ConnectorConstraint, CouplingConstraint, TieConstraint
from femaster_api.model.loads import InertialLoad, NodalForce, PressureLoad, SurfaceTraction, ThermalLoad, VolumeLoad
from femaster_api.model.sections import BeamSection, ShellSection, SolidSection, TrussSection
from femaster_api.model.sets import EntitySet, EntityType
from femaster_api.validation.diagnostics import Diagnostics


def validate_model(model) -> Diagnostics:
    diagnostics = Diagnostics()

    nodes = set(model.nodes.all())
    elements = set(model.elements.all())

    if not nodes:
        diagnostics.error("model has no nodes", location="nodes")
    if not elements:
        diagnostics.warning("model has no elements", location="elements")

    for element in model.elements.all():
        for node in element.nodes:
            if node not in nodes:
                diagnostics.error(f"element {element.id} references a node outside this model", location="elements")

    for item in model.sets.all():
        for member in item.members:
            if item.entity_type == EntityType.NODE and not model.nodes.has(member):
                diagnostics.error(f"set {item.name!r} references a node outside this model", location="sets")
            elif item.entity_type == EntityType.ELEMENT and not model.elements.has(member):
                diagnostics.error(f"set {item.name!r} references an element outside this model", location="sets")

    for section in model.sections.all():
        if section.material and not model.materials.has(section.material):
            diagnostics.error(f"section {section.name!r} references material outside this model", location="sections")
        if not model.sets.has(section.element_set, EntityType.ELEMENT):
            diagnostics.error(f"section {section.name!r} references element set outside this model", location="sections")
        if isinstance(section, BeamSection) and not model.sections.has_profile(section.profile):
            diagnostics.error(f"beam section {section.name!r} references profile outside this model", location="sections")
        if isinstance(section, (SolidSection, ShellSection)) and section.orientation and not model.orientations.has(section.orientation):
            diagnostics.error(f"section {section.name!r} references orientation outside this model", location="sections")
        if isinstance(section, TrussSection) and section.area <= 0:
            diagnostics.error(f"truss section {section.name!r} area must be positive", location="sections")

    for load in model.loads.all():
        _check_load(model, diagnostics, load)

    for support in model.supports.all():
        if support.orientation and not model.orientations.has(support.orientation):
            diagnostics.error("support references orientation outside this model", location="supports")
        _check_target(model, diagnostics, support.target, location="supports")

    for collector in model.support_collectors.all():
        for support in collector.supports:
            if support not in model.supports.all():
                diagnostics.error(f"support collector {collector.name!r} references a support outside this model", location="support_collectors")

    for step in model.steps.all():
        for collector in step.loads:
            if not model.load_collectors.has(collector):
                diagnostics.error(f"step {step.name!r} references load collector outside this model", location="steps")
        for collector in step.supports:
            if not model.support_collectors.has(collector):
                diagnostics.error(f"step {step.name!r} references support collector outside this model", location="steps")
        density = getattr(step, "density", None)
        orientation = getattr(step, "orientation", None)
        initial_velocity = getattr(step, "initial_velocity", None)
        if density and not model.fields.has(density):
            diagnostics.error(f"step {step.name!r} references density field outside this model", location="steps")
        if orientation and not model.fields.has(orientation):
            diagnostics.error(f"step {step.name!r} references orientation field outside this model", location="steps")
        if initial_velocity and not model.fields.has(initial_velocity):
            diagnostics.error(f"step {step.name!r} references initial velocity field outside this model", location="steps")

    for constraint in model.constraints.all():
        if isinstance(constraint, ConnectorConstraint) and not model.orientations.has(constraint.coordinate_system):
            diagnostics.error("connector references coordinate system outside this model", location="constraints")
        if isinstance(constraint, CouplingConstraint):
            if not model.sets.has(constraint.master, EntityType.NODE):
                diagnostics.error("coupling references master node set outside this model", location="constraints")
            if not model.sets.has(constraint.slave, constraint.slave.entity_type):
                diagnostics.error("coupling references slave set outside this model", location="constraints")
        if isinstance(constraint, TieConstraint):
            if not model.sets.has(constraint.master, EntityType.SURFACE):
                diagnostics.warning(f"tie master {constraint.master.name!r} is not a known surface set", location="constraints")

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
        if not model.sets.has(load.target, EntityType.ELEMENT):
            diagnostics.error("inertial load references element set outside this model", location="loads")


def _check_target(model, diagnostics: Diagnostics, target: object, *, location: str) -> None:
    if target is None:
        return
    if isinstance(target, EntitySet) and not model.sets.has(target, target.entity_type):
        diagnostics.error("target references a set outside this model", location=location)
