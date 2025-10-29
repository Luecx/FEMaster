from __future__ import annotations

from typing import Iterable, List, Optional, Union

from ..loads.load_collector import LoadCollector
from ..supports.support_collector import SupportCollector

LoadCollectorInput = Union[LoadCollector, Iterable[LoadCollector]]
SupportCollectorInput = Union[SupportCollector, Iterable[SupportCollector]]


class LinearBucklingStep:
    """Linear buckling analysis step with load/support collectors and eigenvalue count."""

    loadcase_type = "LINEAR BUCKLING"

    def __init__(
        self,
        num_eigenvalues: int,
        sigma: int = None,
        load_collectors: Optional[LoadCollectorInput] = None,
        support_collectors: Optional[SupportCollectorInput] = None,
    ) -> None:
        if num_eigenvalues <= 0:
            raise ValueError("num_eigenvalues must be positive.")
        self.num_eigenvalues = int(num_eigenvalues)
        self.sigma = int(sigma)

        self.load_collectors: List[LoadCollector] = []
        self.support_collectors: List[SupportCollector] = []

        if load_collectors is not None:
            self.extend_load_collectors(load_collectors)
        if support_collectors is not None:
            self.extend_support_collectors(support_collectors)

    def add_load_collector(self, collector: LoadCollector) -> LoadCollector:
        if not isinstance(collector, LoadCollector):
            raise TypeError("LinearBucklingStep only accepts LoadCollector instances.")
        self.load_collectors.append(collector)
        return collector

    def extend_load_collectors(self, collectors: LoadCollectorInput) -> None:
        if isinstance(collectors, LoadCollector):
            self.add_load_collector(collectors)
            return
        for collector in collectors:
            self.add_load_collector(collector)

    def add_support_collector(self, collector: SupportCollector) -> SupportCollector:
        if not isinstance(collector, SupportCollector):
            raise TypeError("LinearBucklingStep only accepts SupportCollector instances.")
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

        blocks.append("*LOADS")
        for collector in self.load_collectors:
            blocks.append(collector.name)

        blocks.append("*SUPPORTS")
        for collector in self.support_collectors:
            blocks.append(collector.name)

        blocks.append("*NUMEIGENVALUES")
        blocks.append(str(self.num_eigenvalues))

        if self.sigma:
            blocks.append("*SIGMA")
            blocks.append(str(self.sigma))

        return "\n".join(blocks)

    def to_asami(self) -> str:
        raise NotImplementedError("LinearBucklingStep.to_asami not implemented yet.")
