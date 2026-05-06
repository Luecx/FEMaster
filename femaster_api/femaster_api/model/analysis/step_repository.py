"""Repository for analysis steps."""

from __future__ import annotations

from typing import Iterator

from .buckling_step import BucklingStep
from .modal_step import ModalStep
from .static_step import StaticStep
from .topology_static_step import TopologyStaticStep
from .transient_step import TransientStep

Step = StaticStep | TopologyStaticStep | ModalStep | BucklingStep | TransientStep


class StepRepository:
    """Repository for FEMaster loadcases."""

    def __init__(self) -> None:
        self._steps: list[Step] = []

    def add(self, step: Step) -> Step:
        if not isinstance(step, (StaticStep, TopologyStaticStep, ModalStep, BucklingStep, TransientStep)):
            raise TypeError("step must be an analysis step object")
        self._steps.append(step)
        return step

    def all(self) -> tuple[Step, ...]:
        return tuple(self._steps)

    def __getitem__(self, index: int | slice) -> Step | tuple[Step, ...]:
        if isinstance(index, slice):
            return tuple(self._steps[index])
        return self._steps[index]

    def __iter__(self) -> Iterator[Step]:
        return iter(self._steps)

    def __len__(self) -> int:
        return len(self._steps)
