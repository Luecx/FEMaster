# fembridge/sets/nodesets.py
from __future__ import annotations
from typing import List, Optional, Union, Iterator
from .nodeset import NodeSet

class NodeSets:
    def __init__(self) -> None:
        self._items: List[Optional[NodeSet]] = []

    def add(self, s: NodeSet) -> NodeSet:
        self._items.append(s)
        return s

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self) -> Iterator[NodeSet]:
        for s in self._items:
            if s is not None:
                yield s

    def __getitem__(self, key: Union[int, str]) -> NodeSet:
        """
        - nodesets[idx] -> NodeSet by index
        - nodesets["NAME"] -> NodeSet by its .name (first match)
        """
        if isinstance(key, int):
            val = self._items[key]
            if val is None:
                raise KeyError(f"NodeSet at index {key} is None.")
            return val
        if isinstance(key, str):
            for s in self._items:
                if s is not None and s.name == key:
                    return s
            raise KeyError(f"NodeSet with name '{key}' not found.")
        raise TypeError("NodeSets.__getitem__ expects int or str.")

    def __contains__(self, key: Union[int, str]) -> bool:
        if isinstance(key, int):
            return 0 <= key < len(self._items) and self._items[key] is not None
        if isinstance(key, str):
            return any(s is not None and s.name == key for s in self._items)
        return False

    def to_femaster(self) -> str:
        blocks = []
        for s in self._items:
            if s is None:
                continue
            blocks.append(s.to_femaster())
        return "\n".join(blocks)

    def to_asami(self) -> str:
        raise NotImplementedError("NodeSets.to_asami not implemented yet.")
