"""General shell ABD elasticity."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .elasticity import Elasticity


class ABDElasticity(Elasticity):
    """General shell ABD constitutive values in native input order."""

    def __init__(self, values: Iterable[float]) -> None:
        self.values = tuple(float(value) for value in values)

    def export(self) -> str:
        lines = [keyword("ELASTIC", TYPE="ABD")]
        for start in range(0, len(self.values), 8):
            lines.append(csv(self.values[start:start + 8]))
        return block(lines)
