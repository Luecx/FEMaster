"""Eigenfrequency analysis step."""

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
            *self.common_export_lines(include_loads=False),
            keyword("NUMEIGENVALUES"),
            csv((self.number_of_modes,)),
        ])
