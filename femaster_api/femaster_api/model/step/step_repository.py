"""Named repository of analysis steps in execution order.

Step names provide semantic access while insertion order defines the order in
which loadcases appear in the exported deck.  Procedure-specific syntax is
entirely delegated to the concrete ``step_*`` classes.
"""

from __future__ import annotations

from ..common.named_repository import NamedRepository
from .step import Step


class StepRepository(NamedRepository[Step]):
    """Own named analysis procedures in execution order."""

    def export(self) -> str:
        return "\n\n".join(step.export() for step in self)
