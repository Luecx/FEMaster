"""Undamped eigenfrequency extraction with concrete support collectors.

A modal step does not activate structural load collectors.  It stores the actual
``SupportCollector`` objects that define its boundary conditions and requests a
fixed number of eigenmodes.  Native collector names are emitted only by ``Step``.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from ..support.support_collector import SupportCollector
from .step import Step


class ModalStep(Step):
    """Undamped eigenfrequency extraction."""

    type_name = "EIGENFREQ"

    def __init__(
        self,
        name: str,
        number_of_modes: int,
        *,
        supports: Iterable[SupportCollector] = (),
    ) -> None:
        super().__init__(name, supports=supports)
        self.number_of_modes = int(number_of_modes)

    def export(self) -> str:
        return block([
            *self._common_lines(include_loads=False),
            keyword("NUMEIGENVALUES"),
            csv((self.number_of_modes,)),
        ])
