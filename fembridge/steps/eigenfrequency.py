from __future__ import annotations

from typing import Iterable, List, Optional, Union

from ..supports.support_collector import SupportCollector

SupportCollectorInput = Union[SupportCollector, Iterable[SupportCollector]]


class EigenfrequencyStep:
    """Eigenfrequency analysis step with support collectors and target eigenvalue count."""

    loadcase_type = "EIGENFREQ"

    def __init__(
        self,
        num_eigenvalues: int,
        support_collectors: Optional[SupportCollectorInput] = None,
    ) -> None:
        if num_eigenvalues <= 0:
            raise ValueError("num_eigenvalues must be positive.")
        self.num_eigenvalues = int(num_eigenvalues)

        self.support_collectors: List[SupportCollector] = []
        if support_collectors is not None:
            self.extend_support_collectors(support_collectors)

    def add_support_collector(self, collector: SupportCollector) -> SupportCollector:
        if not isinstance(collector, SupportCollector):
            raise TypeError("EigenfrequencyStep only accepts SupportCollector instances.")
        self.support_collectors.append(collector)
        return collector

    def extend_support_collectors(self, collectors: SupportCollectorInput) -> None:
        if isinstance(collectors, SupportCollector):
            self.add_support_collector(collectors)
            return
        for collector in collectors:
            self.add_support_collector(collector)

    def to_femaster(self) -> str:
        blocks: List[str] = [f"*LOADCASE, TYPE={self.loadcase_type}"]

        blocks.append("*SUPPORTS")
        for collector in self.support_collectors:
            blocks.append(collector.name)

        blocks.append("*NUMEIGENVALUES")
        blocks.append(str(self.num_eigenvalues))

        return "\n".join(blocks)

    def to_asami(self) -> str:
        raise NotImplementedError("EigenfrequencyStep.to_asami not implemented yet.")
