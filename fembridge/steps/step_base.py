from __future__ import annotations

from typing import Protocol


class StepBase(Protocol):
    """Minimal interface every analysis step should implement."""

    loadcase_type: str

    def to_femaster(self) -> str:
        ...

    def to_asami(self) -> str:
        ...
