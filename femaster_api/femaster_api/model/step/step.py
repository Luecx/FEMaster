"""Base class for one named FEMaster analysis step.

A step corresponds to one native ``*LOADCASE`` and references load/support
collectors by semantic name.  Shared controls such as solver selection and
constraint treatment are emitted here because they are structurally identical
across several concrete analyses.  Procedure-specific commands remain on the
derived ``step_*`` classes.

File naming inside this package intentionally starts with ``step_`` for concrete
procedures so related analysis types stay grouped in directory listings.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import csv, keyword
from ..common.named_object import NamedObject
from .util.constraint_method import ConstraintMethod
from .util.solver_control import SolverControl


class Step(NamedObject):
    """Base class for one named FEMaster loadcase/analysis procedure."""

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
        self.loads = tuple(str(item) for item in loads)
        self.supports = tuple(str(item) for item in supports)
        self.solver = solver
        self.constraint_method = constraint_method

    def _common_lines(self, *, include_loads: bool = True) -> list[str]:
        """Build the common loadcase header and shared child controls."""

        lines = [keyword("LOADCASE", TYPE=self.type_name, NAME=self.name)]

        if self.supports:
            lines.extend([keyword("SUPPORTS"), csv(self.supports)])

        if include_loads and self.loads:
            lines.extend([keyword("LOADS"), csv(self.loads)])

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
