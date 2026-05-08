"""Repository for analysis loadcases."""

from __future__ import annotations

from dataclasses import replace
from typing import Iterator

from femaster_api.model.names import require_name

from .eigenfrequency_loadcase import EigenfrequencyLoadcase
from .linear_buckling_loadcase import LinearBucklingLoadcase
from .linear_static_loadcase import LinearStaticLoadcase
from .linear_transient_loadcase import LinearTransientLoadcase
from .topology_static_loadcase import TopologyStaticLoadcase

Loadcase = (
    LinearStaticLoadcase
    | TopologyStaticLoadcase
    | EigenfrequencyLoadcase
    | LinearBucklingLoadcase
    | LinearTransientLoadcase
)


class LoadcaseRepository:
    """Named container for analysis loadcases."""

    def __init__(self) -> None:
        self._loadcases: dict[str, Loadcase] = {}

    def add(self, loadcase: Loadcase) -> Loadcase:
        """Add a loadcase and return the stored object."""

        if not isinstance(
            loadcase,
            (
                LinearStaticLoadcase,
                TopologyStaticLoadcase,
                EigenfrequencyLoadcase,
                LinearBucklingLoadcase,
                LinearTransientLoadcase,
            ),
        ):
            raise TypeError("loadcase must be an analysis loadcase object")
        item = replace(loadcase, name=require_name(loadcase.name))
        self._loadcases[item.name] = item
        return item

    def all(self) -> tuple[Loadcase, ...]:
        return tuple(self._loadcases[key] for key in sorted(self._loadcases))

    def get(self, name: str) -> Loadcase:
        key = require_name(name)
        try:
            return self._loadcases[key]
        except KeyError as exc:
            raise KeyError(f"unknown loadcase: {key}") from exc

    def __getitem__(self, name: str) -> Loadcase:
        return self.get(name)

    def __contains__(self, value: object) -> bool:
        if isinstance(value, str):
            return require_name(value) in self._loadcases
        return any(item is value for item in self._loadcases.values())

    def __iter__(self) -> Iterator[Loadcase]:
        return iter(self.all())

    def __len__(self) -> int:
        return len(self._loadcases)
