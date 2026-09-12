"""Prescribed structural degrees of freedom on a node target.

A ``Support`` stores up to six translational/rotational values in FEMaster DOF
order.  ``None`` means that a component is unconstrained, while a numeric value
prescribes the corresponding DOF.  The object belongs directly to one
``SupportCollector`` and receives that collector name only during export.

Trailing unconstrained components are omitted from the native data row without
changing the position of interior ``None`` entries.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.typing import EntityReference


class Support:
    """One structural support definition owned by a support collector."""

    def __init__(
        self,
        target: EntityReference,
        values: Iterable[float | None],
        *,
        orientation: str | None = None,
    ) -> None:
        self.target = target
        self.values = tuple(
            None if value is None else float(value)
            for value in values
        )
        self.orientation = orientation

        if len(self.values) > 6:
            raise ValueError("Support accepts at most 6 prescribed DOF values")

    def export(self, collector: str) -> str:
        """Export this support with the owning collector name."""

        values = list(self.values)
        while values and values[-1] is None:
            values.pop()

        return block([
            keyword(
                "SUPPORT",
                SUPPORT_COLLECTOR=collector,
                ORIENTATION=self.orientation,
            ),
            csv((self.target, *values)),
        ])
