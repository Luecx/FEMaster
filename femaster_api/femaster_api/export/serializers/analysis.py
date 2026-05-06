"""FEMaster serialization for analysis steps."""

from __future__ import annotations

from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.analysis import (
    BucklingStep,
    ModalStep,
    StaticStep,
    StepRepository,
    TopologyStaticStep,
    TransientStep,
)


def write_steps(steps: StepRepository) -> str:
    blocks = [_write_step(step) for step in steps.all()]
    return "\n\n".join(block for block in blocks if block)


def _write_step(step) -> str:
    if isinstance(step, StaticStep):
        return _write_static(step)
    if isinstance(step, TopologyStaticStep):
        return _write_topology_static(step)
    if isinstance(step, ModalStep):
        return _write_modal(step)
    if isinstance(step, BucklingStep):
        return _write_buckling(step)
    if isinstance(step, TransientStep):
        return _write_transient(step)
    raise TypeError(f"unsupported analysis step: {type(step).__name__}")


def _write_static(step: StaticStep) -> str:
    lines = [keyword("LOADCASE", TYPE="LINEARSTATIC", NAME=step.name)]
    _write_supports(lines, step.supports)
    _write_loads(lines, step.loads)
    _write_solver(lines, step)
    _write_constraint_method(lines, step)
    if step.inertia_relief:
        lines.append(keyword("INERTIARELIEF", CONSIDER_POINT_MASSES=int(step.inertia_relief_consider_point_masses)))
    if step.rebalance_loads:
        lines.append(keyword("REBALANCELOADS"))
    if step.request_stiffness is not None:
        lines.append(keyword("REQUESTSTIFFNESS", FILE=step.request_stiffness))
    if step.constraint_summary:
        lines.append(keyword("CONSTRAINTSUMMARY"))
    return block(lines)


def _write_topology_static(step: TopologyStaticStep) -> str:
    lines = [keyword("LOADCASE", TYPE="LINEARSTATICTOPO", NAME=step.name)]
    _write_supports(lines, step.supports)
    _write_loads(lines, step.loads)
    _write_solver(lines, step)
    _write_constraint_method(lines, step)
    lines.append(keyword("TOPODENSITY", FIELD=step.density.name))
    if step.orientation is not None:
        lines.append(keyword("TOPOORIENT", FIELD=step.orientation.name))
    if step.exponent is not None:
        lines.append(keyword("TOPOEXPONENT"))
        lines.append(csv((step.exponent,)))
    return block(lines)


def _write_modal(step: ModalStep) -> str:
    lines = [keyword("LOADCASE", TYPE="EIGENFREQ", NAME=step.name)]
    _write_supports(lines, step.supports)
    _write_constraint_method(lines, step)
    lines.append(keyword("NUMEIGENVALUES"))
    lines.append(csv((step.number_of_modes,)))
    return block(lines)


def _write_buckling(step: BucklingStep) -> str:
    lines = [keyword("LOADCASE", TYPE="LINEARBUCKLING", NAME=step.name)]
    _write_supports(lines, step.supports)
    _write_loads(lines, step.loads)
    _write_solver(lines, step)
    _write_constraint_method(lines, step)
    lines.append(keyword("NUMEIGENVALUES"))
    lines.append(csv((step.number_of_modes,)))
    if step.sigma is not None:
        lines.append(keyword("SIGMA"))
        lines.append(csv((step.sigma,)))
    if step.request_stiffness is not None:
        lines.append(keyword("REQUESTSTIFFNESS", FILE=step.request_stiffness))
    if step.request_stgeom is not None:
        lines.append(keyword("REQUESTSTGEOM", FILE=step.request_stgeom))
    return block(lines)


def _write_transient(step: TransientStep) -> str:
    lines = [keyword("LOADCASE", TYPE="LINEARTRANSIENT", NAME=step.name)]
    _write_supports(lines, step.supports)
    _write_loads(lines, step.loads)
    _write_solver(lines, step)
    _write_constraint_method(lines, step)
    lines.append(keyword("TIME"))
    lines.append(csv((step.time.start, step.time.end, step.time.step)))
    if step.newmark is not None:
        lines.append(keyword("NEWMARK"))
        lines.append(csv((step.newmark.beta, step.newmark.gamma)))
    if step.damping is not None:
        lines.append(keyword("DAMPING", TYPE="RAYLEIGH"))
        lines.append(csv((step.damping.alpha, step.damping.beta)))
    if step.write_every is not None:
        lines.append(keyword("WRITEEVERY", TYPE=step.write_every_type or "STEPS"))
        lines.append(csv((step.write_every,)))
    if step.initial_velocity is not None:
        lines.append(keyword("INITIALVELOCITY", FIELD=step.initial_velocity.name))
    return block(lines)


def _write_supports(lines: list[str], supports) -> None:
    if supports:
        lines.append(keyword("SUPPORTS"))
        lines.append(csv(collector.name for collector in supports))


def _write_loads(lines: list[str], loads) -> None:
    if loads:
        lines.append(keyword("LOADS"))
        lines.append(csv(collector.name for collector in loads))


def _write_solver(lines: list[str], step) -> None:
    if step.solver is not None:
        lines.append(keyword("SOLVER", DEVICE=step.solver.device.value, METHOD=step.solver.method.value))


def _write_constraint_method(lines: list[str], step) -> None:
    if step.constraint_method is not None:
        lines.append(keyword("CONSTRAINTMETHOD", TYPE=step.constraint_method.value))
