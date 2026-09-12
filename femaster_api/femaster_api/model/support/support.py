"""Prescribed structural support."""

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..common.typing import EntityReference


class Support:
    """Prescribed translational/rotational DOFs on a node or node region."""

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
