"""Base analysis step."""

from collections.abc import Iterable

from ..common.format import csv, keyword
from ..common.named_object import NamedObject
from .constraint_method import ConstraintMethod
from .solver_control import SolverControl


class Step(NamedObject):
    """Base named FEMaster LOADCASE."""

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

    def common_export_lines(self, *, include_loads: bool = True) -> list[str]:
        lines = [keyword("LOADCASE", TYPE=self.type_name, NAME=self.name)]
        if self.supports:
            lines.extend([keyword("SUPPORTS"), csv(self.supports)])
        if include_loads and self.loads:
            lines.extend([keyword("LOADS"), csv(self.loads)])
        if self.solver is not None:
            lines.append(self.solver.export())
        if self.constraint_method is not None:
            lines.append(keyword("CONSTRAINTMETHOD", TYPE=self.constraint_method.value))
        return lines

    def export(self) -> str:
        raise NotImplementedError
