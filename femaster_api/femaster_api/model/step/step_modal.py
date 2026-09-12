"""Undamped eigenfrequency extraction step.

A modal step does not activate structural load collectors; it uses support
collectors and requests a fixed number of eigenmodes.  The concrete class owns
the ``*NUMEIGENVALUES`` child block while shared loadcase syntax remains in
``Step``.
"""

from __future__ import annotations

from collections.abc import Iterable

from ..common.format import block, csv, keyword
from .step import Step


class ModalStep(Step):
    """Undamped eigenfrequency extraction."""

    type_name = "EIGENFREQ"

    def __init__(
        self,
        name: str,
        number_of_modes: int,
        *,
        supports: Iterable[str] = (),
    ) -> None:
        super().__init__(name, supports=supports)
        self.number_of_modes = int(number_of_modes)

    def export(self) -> str:
        return block([
            *self._common_lines(include_loads=False),
            keyword("NUMEIGENVALUES"),
            csv((self.number_of_modes,)),
        ])
