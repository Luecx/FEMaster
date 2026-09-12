"""Base class for one named FEMaster analysis step with object references.

A step corresponds to one native ``*LOADCASE`` and stores concrete
``LoadCollector`` / ``SupportCollector`` objects.  Collector names exist only in
the exported input syntax; callers cannot pass strings in place of the objects
that the analysis actually activates.

Shared numerical controls are emitted here because they are structurally common
across several procedures.  Procedure-specific controls remain on the concrete
``step_*`` classes.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import csv, keyword
from ..common.named_object import NamedObject
from ..load.load_collector import LoadCollector
from ..support.support_collector import SupportCollector
from .util.constraint_method import ConstraintMethod
from .util.solver_control import SolverControl


class Step(NamedObject):
    """Base class for one named FEMaster loadcase/analysis procedure."""

    type_name = "LOADCASE"

    def __init__(
        self,
        name: str,
        *,
        loads: Iterable[LoadCollector] = (),
        supports: Iterable[SupportCollector] = (),
        solver: SolverControl | None = None,
        constraint_method: ConstraintMethod | None = None,
    ) -> None:
        super().__init__(name)
        self.loads = tuple(loads)
        self.supports = tuple(supports)

        if not all(isinstance(item, LoadCollector) for item in self.loads):
            raise TypeError("loads must contain LoadCollector objects")
        if not all(isinstance(item, SupportCollector) for item in self.supports):
            raise TypeError("supports must contain SupportCollector objects")
        if solver is not None and not isinstance(solver, SolverControl):
            raise TypeError("solver must be a SolverControl object")
        if constraint_method is not None and not isinstance(
            constraint_method,
            ConstraintMethod,
        ):
            raise TypeError("constraint_method must be ConstraintMethod")

        self.solver = solver
        self.constraint_method = constraint_method

    def _common_lines(self, *, include_loads: bool = True) -> list[str]:
        """Build the common loadcase header and shared child controls."""

        lines = [keyword("LOADCASE", TYPE=self.type_name, NAME=self.name)]

        if self.supports:
            lines.extend([
                keyword("SUPPORTS"),
                csv(tuple(item.name for item in self.supports)),
            ])

        if include_loads and self.loads:
            lines.extend([
                keyword("LOADS"),
                csv(tuple(item.name for item in self.loads)),
            ])

        if self.solver is not None:
            lines.append(self.solver.export())

        if self.constraint_method is not None:
            lines.append(
                keyword(
                    "CONSTRAINTMETHOD",
                    TYPE=self.constraint_method.value,
                )
            )

        return lines

    def export(self) -> str:
        """Export the complete procedure-specific loadcase block."""

        raise NotImplementedError
