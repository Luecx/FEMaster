from __future__ import annotations

from typing import Iterable, Iterator, List

from .step_base import StepBase


class Steps:
    """Container managing the analysis steps belonging to a model."""

    def __init__(self, entries: Iterable[StepBase] = ()) -> None:
        self._items: List[StepBase] = []
        self.extend(entries)

    def add(self, step: StepBase) -> StepBase:
        if not hasattr(step, "to_femaster"):
            raise TypeError("Steps can only contain objects providing to_femaster().")
        self._items.append(step)
        return step

    def extend(self, steps: Iterable[StepBase]) -> None:
        for step in steps:
            self.add(step)

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self) -> Iterator[StepBase]:
        return iter(self._items)

    def __getitem__(self, index: int) -> StepBase:
        return self._items[index]

    def to_femaster(self) -> str:
        blocks = [step.to_femaster() for step in self._items]
        return "\n".join(block for block in blocks if block)

    def to_asami(self) -> str:
        raise NotImplementedError("Steps.to_asami not implemented yet.")
