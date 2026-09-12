"""Repository of analysis steps."""

from ..common.named_repository import NamedRepository
from .step import Step


class StepRepository(NamedRepository[Step]):
    """Named analysis-step repository preserving execution order."""

    def export(self) -> str:
        return "\n\n".join(step.export() for step in self)
