"""Named FEMaster analysis-step definitions.

Every concrete step owns the keyword sequence needed for its analysis type.
Shared controls remain small value objects and are rendered by the Step base
class. StepRepository only provides ordered/name-based ownership and export.
"""

from __future__ import annotations

from enum import Enum
from typing import Iterable

from .._format import block, csv, keyword
from ..repository import NamedObject, NamedRepository


class SolverDevice(Enum):
    CPU = "CPU"
    GPU = "GPU"


class SolverMethod(Enum):
    DIRECT   = "DIRECT"
    INDIRECT = "INDIRECT"


class ConstraintMethod(Enum):
    NULLSPACE = "NULLSPACE"
    LAGRANGE  = "LAGRANGE"


class SolverControl:
    """Linear solver device/method selection shared by analysis steps."""

    def __init__(
        self,
        device: SolverDevice = SolverDevice.CPU,
        method: SolverMethod = SolverMethod.DIRECT,
    ) -> None:
        self.device = device
        self.method = method

    def to_femaster(self) -> str:
        return keyword("SOLVER", DEVICE=self.device.value, METHOD=self.method.value)


class Step(NamedObject):
    """Base class for one named FEMaster LOADCASE."""

    type_name = "LOADCASE"

    def __init__(
        self,
        name: str,
        *,
        loads: Iterable[str] = (),
        supports: Iterable[str] = (),
        solver: SolverControl | None = None,
        constraint_method: ConstraintMethod | None = None,
    ) -> None:
        super().__init__(name)
        self.loads             = tuple(str(item) for item in loads)
        self.supports          = tuple(str(item) for item in supports)
        self.solver            = solver
        self.constraint_method = constraint_method

    def _common_lines(self, *, include_loads: bool = True) -> list[str]:
        """Build common collector and solver controls in FEMaster order."""

        lines = [keyword("LOADCASE", TYPE=self.type_name, NAME=self.name)]

        if self.supports:
            lines.extend([keyword("SUPPORTS"), csv(self.supports)])

        if include_loads and self.loads:
            lines.extend([keyword("LOADS"), csv(self.loads)])

        if self.solver is not None:
            lines.append(self.solver.to_femaster())

        if self.constraint_method is not None:
            lines.append(keyword("CONSTRAINTMETHOD", TYPE=self.constraint_method.value))

        return lines

    def to_femaster(self) -> str:
        raise NotImplementedError


class StaticStep(Step):
    """Linear static structural analysis."""

    type_name = "LINEARSTATIC"

    def __init__(
        self,
        name: str,
        *,
        loads: Iterable[str] = (),
        supports: Iterable[str] = (),
        solver: SolverControl | None = None,
        constraint_method: ConstraintMethod | None = None,
        inertia_relief: bool = False,
        rebalance_loads: bool = False,
    ) -> None:
        super().__init__(
            name,
            loads=loads,
            supports=supports,
            solver=solver,
            constraint_method=constraint_method,
        )
        self.inertia_relief  = bool(inertia_relief)
        self.rebalance_loads = bool(rebalance_loads)

    def to_femaster(self) -> str:
        lines = self._common_lines()
        if self.inertia_relief:
            lines.append(keyword("INERTIARELIEF"))
        if self.rebalance_loads:
            lines.append(keyword("REBALANCELOADS"))
        return block(lines)


class ModalStep(Step):
    """Undamped eigenfrequency extraction."""

    type_name = "EIGENFREQ"

    def __init__(self, name: str, number_of_modes: int, *, supports: Iterable[str] = ()) -> None:
        super().__init__(name, supports=supports)
        self.number_of_modes = int(number_of_modes)

    def to_femaster(self) -> str:
        return block([
            *self._common_lines(include_loads=False),
            keyword("NUMEIGENVALUES"),
            csv((self.number_of_modes,)),
        ])


class BucklingStep(Step):
    """Linearized eigenvalue buckling analysis around a static preload."""

    type_name = "LINEARBUCKLING"

    def __init__(
        self,
        name: str,
        number_of_modes: int,
        *,
        loads: Iterable[str] = (),
        supports: Iterable[str] = (),
        solver: SolverControl | None = None,
        sigma: float | None = None,
    ) -> None:
        super().__init__(name, loads=loads, supports=supports, solver=solver)
        self.number_of_modes = int(number_of_modes)
        self.sigma           = None if sigma is None else float(sigma)

    def to_femaster(self) -> str:
        lines = [
            *self._common_lines(),
            keyword("NUMEIGENVALUES"),
            csv((self.number_of_modes,)),
        ]
        if self.sigma is not None:
            lines.extend([keyword("SIGMA"), csv((self.sigma,))])
        return block(lines)


class NonlinearStaticStep(Step):
    """Geometrically/materially nonlinear static analysis."""

    type_name = "NONLINEARSTATIC"

    def __init__(
        self,
        name: str,
        *,
        loads: Iterable[str] = (),
        supports: Iterable[str] = (),
        solver: SolverControl | None = None,
        constraint_method: ConstraintMethod | None = None,
        control: str = "LOAD",
        increments: int | None = None,
        max_increments: int | None = None,
        initial_increment: float | None = None,
        minimum_increment: float | None = None,
        maximum_increment: float | None = None,
        max_iterations: int | None = None,
        tolerance: float | None = None,
    ) -> None:
        super().__init__(
            name,
            loads=loads,
            supports=supports,
            solver=solver,
            constraint_method=constraint_method,
        )
        self.control           = str(control).upper()
        self.increments        = increments
        self.max_increments    = max_increments
        self.initial_increment = initial_increment
        self.minimum_increment = minimum_increment
        self.maximum_increment = maximum_increment
        self.max_iterations    = max_iterations
        self.tolerance         = tolerance

    def to_femaster(self) -> str:
        keys = {
            "CONTROL": self.control,
            "INCREMENTS": self.increments,
            "MAX_INCREMENTS": self.max_increments,
            "INITIAL_INCREMENT": self.initial_increment,
            "MINIMUM_INCREMENT": self.minimum_increment,
            "MAXIMUM_INCREMENT": self.maximum_increment,
            "MAXITER": self.max_iterations,
            "TOL": self.tolerance,
        }
        return block([
            *self._common_lines(),
            keyword("NONLINEAR", **keys),
        ])


class TimeControl:
    """Start/end/step definition for transient integration."""

    def __init__(self, start: float, end: float, step: float) -> None:
        self.start = float(start)
        self.end   = float(end)
        self.step  = float(step)


class NewmarkControl:
    """Newmark beta/gamma integration parameters."""

    def __init__(self, beta: float = 0.25, gamma: float = 0.5) -> None:
        self.beta  = float(beta)
        self.gamma = float(gamma)


class RayleighDamping:
    """Mass- and stiffness-proportional Rayleigh damping coefficients."""

    def __init__(self, alpha: float = 0.0, beta: float = 0.0) -> None:
        self.alpha = float(alpha)
        self.beta  = float(beta)


class TransientStep(Step):
    """Linear structural transient analysis with Newmark integration."""

    type_name = "LINEARTRANSIENT"

    def __init__(
        self,
        name: str,
        time: TimeControl,
        *,
        loads: Iterable[str] = (),
        supports: Iterable[str] = (),
        solver: SolverControl | None = None,
        newmark: NewmarkControl | None = None,
        damping: RayleighDamping | None = None,
        write_every: int | None = None,
    ) -> None:
        super().__init__(name, loads=loads, supports=supports, solver=solver)
        self.time        = time
        self.newmark     = newmark
        self.damping     = damping
        self.write_every = write_every

    def to_femaster(self) -> str:
        lines = [
            *self._common_lines(),
            keyword("TIME"),
            csv((self.time.start, self.time.end, self.time.step)),
        ]
        if self.newmark is not None:
            lines.extend([keyword("NEWMARK"), csv((self.newmark.beta, self.newmark.gamma))])
        if self.damping is not None:
            lines.extend([
                keyword("DAMPING", TYPE="RAYLEIGH"),
                csv((self.damping.alpha, self.damping.beta)),
            ])
        if self.write_every is not None:
            lines.extend([keyword("WRITEEVERY", TYPE="STEPS"), csv((self.write_every,))])
        return block(lines)


class StepRepository(NamedRepository[Step]):
    """Repository of named analysis steps in execution order."""

    def to_femaster(self) -> str:
        return "\n\n".join(step.to_femaster() for step in self)
