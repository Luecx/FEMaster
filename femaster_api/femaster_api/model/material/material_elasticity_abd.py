"""Generalized shell ABD elasticity stored in native FEMaster input order.

This material form retains the complete sequence supplied by the caller and
writes it over continuation rows below ``*ELASTIC, TYPE=ABD``.  Validation of a
section-level 6x6 ABD plus shear representation belongs to ``ABDShellSection``;
this class only represents the material keyword form.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .material_elasticity import Elasticity


class ABDElasticity(Elasticity):
    """Generalized shell elasticity represented by an ordered coefficient list."""

    def __init__(self, values: Iterable[float]) -> None:
        self.values = tuple(float(value) for value in values)

    def export(self) -> str:
        lines = [keyword("ELASTIC", TYPE="ABD")]
        for start in range(0, len(self.values), 8):
            lines.append(csv(self.values[start:start + 8]))
        return block(lines)
