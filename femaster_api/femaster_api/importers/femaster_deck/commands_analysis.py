"""Analysis/loadcase commands."""

from __future__ import annotations

from femaster_api.model import (
    BucklingStep,
    ConstraintMethod,
    ModalStep,
    NewmarkControl,
    NonlinearStaticStep,
    RayleighDamping,
    SolverControl,
    SolverDevice,
    SolverMethod,
    StaticStep,
    TimeControl,
    TopologyStaticStep,
    TransientStep,
)

from .utils import FEMasterInputError, Header, falsey, normalize_loadcase_type, numbers, parse_csv, parse_header, truthy


LOADCASE_SUBCOMMANDS = {
    "LOADS",
    "SUPPORTS",
    "SOLVER",
    "CONSTRAINTMETHOD",
    "INERTIARELIEF",
    "REBALANCELOADS",
    "REQUESTSTIFFNESS",
    "REQUESTSTGEOM",
    "CONSTRAINTSUMMARY",
    "NUMEIGENVALUES",
    "SIGMA",
    "TOPODENSITY",
    "TOPOORIENT",
    "TOPOORIENTATION",
    "TOPOEXPONENT",
    "NONLINEAR",
    "TIME",
    "NEWMARK",
    "DAMPING",
    "WRITEEVERY",
    "INITIALVELOCITY",
}


def cmd_loadcase(parser, header: Header) -> None:
    lc_type = normalize_loadcase_type(header.params.get("TYPE"))
    state = {
        "name": header.params.get("NAME") or lc_type,
        "loads": (),
        "supports": (),
        "solver": None,
        "constraint_method": ConstraintMethod.NULLSPACE,
        "inertia_relief": False,
        "inertia_relief_consider_point_masses": True,
        "rebalance_loads": False,
        "request_stiffness": None,
        "request_stgeom": None,
        "constraint_summary": False,
        "number_of_modes": 10,
        "sigma": None,
        "density": None,
        "orientation": None,
        "exponent": None,
        "time": None,
        "newmark": None,
        "damping": None,
        "write_every_type": None,
        "write_every": None,
        "initial_velocity": None,
        "nonlinear": {},
    }

    while True:
        line = parser.stream.peek()
        if line is None:
            break
        if not line.startswith("*"):
            raise FEMasterInputError(f"unexpected data inside *LOADCASE: {line!r}")
        next_header = parse_header(line)
        if next_header.keyword.startswith("END"):
            parser.stream.pop()
            break
        if next_header.keyword not in LOADCASE_SUBCOMMANDS:
            if parser.next_header_is_top_level(LOADCASE_SUBCOMMANDS):
                break
            raise FEMasterInputError(f"unsupported LOADCASE command '*{next_header.keyword}'")
        sub = parse_header(parser.stream.pop() or "")
        _parse_loadcase_subcommand(parser, lc_type, state, sub)

    parser.model.steps.add(_make_step(lc_type, state))


def _parse_loadcase_subcommand(parser, lc_type: str, state: dict, header: Header) -> None:
    keyword = header.keyword
    if keyword == "LOADS":
        state["loads"] = tuple(parser.model.load_collectors.get(token) for token in _tokens(parser))
    elif keyword == "SUPPORTS":
        state["supports"] = tuple(parser.model.support_collectors.get(token) for token in _tokens(parser))
    elif keyword == "SOLVER":
        _reject_data(parser, "*SOLVER")
        state["solver"] = SolverControl(
            SolverDevice[header.params.get("DEVICE", "CPU").upper()],
            SolverMethod[header.params.get("METHOD", "DIRECT").upper()],
        )
    elif keyword == "CONSTRAINTMETHOD":
        if lc_type not in {"LINEARSTATIC", "LINEARSTATICTOPO", "NONLINEARSTATIC"}:
            raise FEMasterInputError("*CONSTRAINTMETHOD is not supported for this loadcase type")
        _reject_data(parser, "*CONSTRAINTMETHOD")
        state["constraint_method"] = ConstraintMethod[header.params.get("TYPE", "NULLSPACE").upper()]
    elif keyword == "INERTIARELIEF":
        _require(lc_type == "LINEARSTATIC", "*INERTIARELIEF only valid for LINEARSTATIC")
        _reject_data(parser, "*INERTIARELIEF")
        state["inertia_relief"] = True
        state["inertia_relief_consider_point_masses"] = not falsey(header.params.get("CONSIDER_POINT_MASSES"))
    elif keyword == "REBALANCELOADS":
        _require(lc_type == "LINEARSTATIC", "*REBALANCELOADS only valid for LINEARSTATIC")
        _reject_data(parser, "*REBALANCELOADS")
        state["rebalance_loads"] = True
    elif keyword == "REQUESTSTIFFNESS":
        _reject_data(parser, "*REQUESTSTIFFNESS")
        state["request_stiffness"] = header.params.get("FILE") or "stiffness.txt"
    elif keyword == "REQUESTSTGEOM":
        _reject_data(parser, "*REQUESTSTGEOM")
        state["request_stgeom"] = header.params.get("FILE") or "stgeom.txt"
    elif keyword == "CONSTRAINTSUMMARY":
        _reject_data(parser, "*CONSTRAINTSUMMARY")
        state["constraint_summary"] = True
    elif keyword == "NUMEIGENVALUES":
        state["number_of_modes"] = int(_one_number(parser, "*NUMEIGENVALUES"))
    elif keyword == "SIGMA":
        state["sigma"] = _one_number(parser, "*SIGMA")
    elif keyword == "TOPODENSITY":
        field_name = header.params.get("FIELD")
        if not field_name:
            raise FEMasterInputError("*TOPODENSITY requires FIELD")
        _reject_data(parser, "*TOPODENSITY")
        state["density"] = parser.model.fields.get(field_name)
    elif keyword in {"TOPOORIENT", "TOPOORIENTATION"}:
        field_name = header.params.get("FIELD")
        if not field_name:
            raise FEMasterInputError(f"*{keyword} requires FIELD")
        _reject_data(parser, f"*{keyword}")
        state["orientation"] = parser.model.fields.get(field_name)
    elif keyword == "TOPOEXPONENT":
        state["exponent"] = _one_number(parser, "*TOPOEXPONENT")
    elif keyword == "NONLINEAR":
        _reject_data(parser, "*NONLINEAR")
        state["nonlinear"] = _nonlinear_values(header)
    elif keyword == "TIME":
        values = [value for line in parser.consume_data_lines() for value in numbers(line)]
        if len(values) == 2:
            state["time"] = TimeControl(0.0, values[0], values[1])
        elif len(values) >= 3:
            state["time"] = TimeControl(values[0], values[1], values[2])
        else:
            raise FEMasterInputError("*TIME requires t_end,dt or t_start,t_end,dt")
    elif keyword == "NEWMARK":
        values = [value for line in parser.consume_data_lines() for value in numbers(line)]
        if len(values) < 2:
            raise FEMasterInputError("*NEWMARK requires beta,gamma")
        state["newmark"] = NewmarkControl(values[0], values[1])
    elif keyword == "DAMPING":
        values = [value for line in parser.consume_data_lines() for value in numbers(line)]
        if len(values) < 2:
            raise FEMasterInputError("*DAMPING requires alpha,beta")
        state["damping"] = RayleighDamping(values[0], values[1])
    elif keyword == "WRITEEVERY":
        state["write_every_type"] = header.params.get("TYPE", "STEPS").upper()
        state["write_every"] = _one_number(parser, "*WRITEEVERY")
    elif keyword == "INITIALVELOCITY":
        field_name = header.params.get("FIELD")
        if not field_name:
            raise FEMasterInputError("*INITIALVELOCITY requires FIELD")
        _reject_data(parser, "*INITIALVELOCITY")
        state["initial_velocity"] = parser.model.fields.get(field_name)


def _make_step(lc_type: str, state: dict):
    if lc_type == "LINEARSTATIC":
        return StaticStep(
            state["name"],
            state["loads"],
            state["supports"],
            state["solver"],
            state["constraint_method"],
            state["inertia_relief"],
            state["inertia_relief_consider_point_masses"],
            state["rebalance_loads"],
            state["request_stiffness"],
            state["constraint_summary"],
        )
    if lc_type == "LINEARSTATICTOPO":
        _require(state["density"] is not None, "LINEARSTATICTOPO requires *TOPODENSITY")
        return TopologyStaticStep(
            state["name"],
            state["density"],
            state["orientation"],
            state["exponent"],
            state["loads"],
            state["supports"],
            state["solver"],
            state["constraint_method"],
        )
    if lc_type == "EIGENFREQ":
        return ModalStep(state["name"], state["supports"], state["number_of_modes"])
    if lc_type == "LINEARBUCKLING":
        return BucklingStep(
            state["name"],
            state["loads"],
            state["supports"],
            state["number_of_modes"],
            state["sigma"],
            state["solver"],
            None,
            state["request_stiffness"],
            state["request_stgeom"],
        )
    if lc_type == "LINEARTRANSIENT":
        _require(state["time"] is not None, "LINEARTRANSIENT requires *TIME")
        return TransientStep(
            state["name"],
            state["loads"],
            state["time"],
            state["supports"],
            state["solver"],
            None,
            state["newmark"],
            state["damping"],
            state["write_every_type"],
            state["write_every"],
            state["initial_velocity"],
        )
    if lc_type == "NONLINEARSTATIC":
        return NonlinearStaticStep(
            state["name"],
            state["loads"],
            state["supports"],
            state["solver"],
            state["constraint_method"],
            request_stiffness=state["request_stiffness"],
            constraint_summary=state["constraint_summary"],
            **state["nonlinear"],
        )
    raise FEMasterInputError(f"unsupported LOADCASE TYPE={lc_type!r}")


def _nonlinear_values(header: Header) -> dict:
    values = {
        "control": header.params.get("CONTROL", "LOAD").upper(),
        "increments": _int_param(header, "INCREMENTS"),
        "max_increments": _int_param(header, "MAX_INCREMENTS"),
        "initial_increment": _float_param(header, "INITIAL_INCREMENT"),
        "minimum_increment": _float_param(header, "MINIMUM_INCREMENT"),
        "maximum_increment": _float_param(header, "MAXIMUM_INCREMENT"),
        "arc_length_psi": _float_param(header, "ARC_LENGTH_PSI"),
        "adaptive": _bool_param(header, "ADAPTIVE"),
        "growth_factor": _float_param(header, "GROWTH_FACTOR"),
        "cutback_factor": _float_param(header, "CUTBACK_FACTOR"),
        "fast_iterations": _int_param(header, "FAST_ITERATIONS"),
        "slow_iterations": _int_param(header, "SLOW_ITERATIONS"),
        "maximum_cutbacks": _int_param(header, "MAXIMUM_CUTBACKS"),
        "max_iterations": _int_param(header, "MAXITER"),
        "tolerance": _float_param(header, "TOL"),
        "regularize_zero_rows": _bool_param(header, "REGULARIZE_ZERO_ROWS"),
        "regularization_alpha": _float_param(header, "REGULARIZATION_ALPHA"),
    }
    return {key: value for key, value in values.items() if value is not None}


def _tokens(parser) -> list[str]:
    return [token for line in parser.consume_data_lines() for token in parse_csv(line) if token]


def _one_number(parser, command: str) -> float:
    values = [value for line in parser.consume_data_lines() for value in numbers(line)]
    if not values:
        raise FEMasterInputError(f"{command} requires one numeric value")
    return values[0]


def _reject_data(parser, command: str) -> None:
    for line in parser.consume_data_lines():
        if line.strip():
            raise FEMasterInputError(f"{command} does not take data lines")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise FEMasterInputError(message)


def _int_param(header: Header, key: str) -> int | None:
    return None if key not in header.params else int(float(header.params[key]))


def _float_param(header: Header, key: str) -> float | None:
    return None if key not in header.params else float(header.params[key])


def _bool_param(header: Header, key: str) -> bool | None:
    if key not in header.params:
        return None
    value = header.params[key]
    if truthy(value):
        return True
    if falsey(value):
        return False
    raise FEMasterInputError(f"{key} must be ON/OFF or TRUE/FALSE")


ANALYSIS_COMMANDS = {
    "LOADCASE": cmd_loadcase,
}
