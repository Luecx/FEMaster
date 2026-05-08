"""FEMaster serialization for analysis loadcases."""

from __future__ import annotations

from femaster_api.export.femaster_format import block, csv, keyword
from femaster_api.model.analysis import (
    EigenfrequencyLoadcase,
    LinearBucklingLoadcase,
    LinearStaticLoadcase,
    LinearTransientLoadcase,
    LoadcaseRepository,
    TopologyStaticLoadcase,
)


def write_loadcases(loadcases: LoadcaseRepository) -> str:
    """Serialize analysis loadcases."""

    blocks = [_write_loadcase(loadcase) for loadcase in loadcases]
    return "\n\n".join(block for block in blocks if block)


def _write_loadcase(loadcase) -> str:
    if isinstance(loadcase, LinearStaticLoadcase):
        return _write_static(loadcase)
    if isinstance(loadcase, TopologyStaticLoadcase):
        return _write_topology_static(loadcase)
    if isinstance(loadcase, EigenfrequencyLoadcase):
        return _write_modal(loadcase)
    if isinstance(loadcase, LinearBucklingLoadcase):
        return _write_buckling(loadcase)
    if isinstance(loadcase, LinearTransientLoadcase):
        return _write_transient(loadcase)
    raise TypeError(f"unsupported analysis loadcase: {type(loadcase).__name__}")


def _write_static(loadcase: LinearStaticLoadcase) -> str:
    lines = [keyword("LOADCASE", TYPE="LINEARSTATIC", NAME=loadcase.name)]
    _write_supports(lines, loadcase.supports)
    _write_loads(lines, loadcase.loads)
    _write_solver(lines, loadcase)
    _write_constraint_method(lines, loadcase)
    if loadcase.inertia_relief:
        lines.append(keyword("INERTIARELIEF", CONSIDER_POINT_MASSES=int(loadcase.inertia_relief_consider_point_masses)))
    if loadcase.rebalance_loads:
        lines.append(keyword("REBALANCELOADS"))
    if loadcase.request_stiffness is not None:
        lines.append(keyword("REQUESTSTIFFNESS", FILE=loadcase.request_stiffness))
    if loadcase.constraint_summary:
        lines.append(keyword("CONSTRAINTSUMMARY"))
    return block(lines)


def _write_topology_static(loadcase: TopologyStaticLoadcase) -> str:
    lines = [keyword("LOADCASE", TYPE="LINEARSTATICTOPO", NAME=loadcase.name)]
    _write_supports(lines, loadcase.supports)
    _write_loads(lines, loadcase.loads)
    _write_solver(lines, loadcase)
    _write_constraint_method(lines, loadcase)
    lines.append(keyword("TOPODENSITY", FIELD=loadcase.density.name))
    if loadcase.orientation is not None:
        lines.append(keyword("TOPOORIENT", FIELD=loadcase.orientation.name))
    if loadcase.exponent is not None:
        lines.append(keyword("TOPOEXPONENT"))
        lines.append(csv((loadcase.exponent,)))
    return block(lines)


def _write_modal(loadcase: EigenfrequencyLoadcase) -> str:
    lines = [keyword("LOADCASE", TYPE="EIGENFREQ", NAME=loadcase.name)]
    _write_supports(lines, loadcase.supports)
    _write_constraint_method(lines, loadcase)
    lines.append(keyword("NUMEIGENVALUES"))
    lines.append(csv((loadcase.number_of_modes,)))
    return block(lines)


def _write_buckling(loadcase: LinearBucklingLoadcase) -> str:
    lines = [keyword("LOADCASE", TYPE="LINEARBUCKLING", NAME=loadcase.name)]
    _write_supports(lines, loadcase.supports)
    _write_loads(lines, loadcase.loads)
    _write_solver(lines, loadcase)
    _write_constraint_method(lines, loadcase)
    lines.append(keyword("NUMEIGENVALUES"))
    lines.append(csv((loadcase.number_of_modes,)))
    if loadcase.sigma is not None:
        lines.append(keyword("SIGMA"))
        lines.append(csv((loadcase.sigma,)))
    if loadcase.request_stiffness is not None:
        lines.append(keyword("REQUESTSTIFFNESS", FILE=loadcase.request_stiffness))
    if loadcase.request_stgeom is not None:
        lines.append(keyword("REQUESTSTGEOM", FILE=loadcase.request_stgeom))
    return block(lines)


def _write_transient(loadcase: LinearTransientLoadcase) -> str:
    lines = [keyword("LOADCASE", TYPE="LINEARTRANSIENT", NAME=loadcase.name)]
    _write_supports(lines, loadcase.supports)
    _write_loads(lines, loadcase.loads)
    _write_solver(lines, loadcase)
    _write_constraint_method(lines, loadcase)
    lines.append(keyword("TIME"))
    lines.append(csv((loadcase.time.start, loadcase.time.end, loadcase.time.step)))
    if loadcase.newmark is not None:
        lines.append(keyword("NEWMARK"))
        lines.append(csv((loadcase.newmark.beta, loadcase.newmark.gamma)))
    if loadcase.damping is not None:
        lines.append(keyword("DAMPING", TYPE="RAYLEIGH"))
        lines.append(csv((loadcase.damping.alpha, loadcase.damping.beta)))
    if loadcase.write_every is not None:
        lines.append(keyword("WRITEEVERY", TYPE=loadcase.write_every_type or "STEPS"))
        lines.append(csv((loadcase.write_every,)))
    if loadcase.initial_velocity is not None:
        lines.append(keyword("INITIALVELOCITY", FIELD=loadcase.initial_velocity.name))
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
