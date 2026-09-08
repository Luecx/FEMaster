"""Lazy mesh extension for the FEMR result reader."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import struct

from .femr_reader import FemrResults, _CHUNK, _DataRef, _decompress, _verify


@dataclass(frozen=True, slots=True)
class FemrElement:
    id: int
    type: str
    node_ids: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class FemrMesh:
    nodes: dict[int, tuple[float, float, float]]
    elements: dict[int, FemrElement]

    def node(self, id: int) -> tuple[float, float, float]:
        return self.nodes[id]

    def element(self, id: int) -> FemrElement:
        return self.elements[id]


class MeshFemrResults(FemrResults):
    """FEMR results whose mesh chunk is loaded only on first access."""

    def __init__(self, path: str | Path) -> None:
        super().__init__(path)
        self._mesh_ref = self._find_mesh()
        self._mesh: FemrMesh | None = None

    @property
    def mesh(self) -> FemrMesh:
        if self._mesh is None:
            if self._mesh_ref is None:
                raise ValueError("FEMR file contains no MESH chunk")
            self._mesh = self._load_mesh(self._mesh_ref)
        return self._mesh

    def _find_mesh(self) -> _DataRef | None:
        self._file.seek(0)
        mesh_ref: _DataRef | None = None
        while header := self._file.read(_CHUNK.size):
            if len(header) != _CHUNK.size:
                raise ValueError("truncated FEMR chunk header")
            kind, compression, stored_size, raw_size, checksum, _ = _CHUNK.unpack(header)
            payload_offset = self._file.tell()
            if kind == b"MESH":
                if mesh_ref is not None:
                    raise ValueError("FEMR file contains multiple MESH chunks")
                mesh_ref = _DataRef(payload_offset, stored_size, raw_size, compression, checksum)
            self._file.seek(stored_size, 1)
        return mesh_ref

    def _load_mesh(self, ref: _DataRef) -> FemrMesh:
        if self._file.closed:
            raise ValueError("cannot load MESH after FEMR file was closed")
        self._file.seek(ref.offset)
        stored = self._file.read(ref.stored_size)
        if len(stored) != ref.stored_size:
            raise ValueError("truncated FEMR MESH payload")
        raw = _decompress(stored, ref.compression, ref.raw_size)
        _verify(raw, ref.checksum, b"MESH")

        if len(raw) < 16:
            raise ValueError("invalid FEMR mesh chunk")
        node_count, element_count = struct.unpack_from("<QQ", raw)
        pos = 16
        nodes: dict[int, tuple[float, float, float]] = {}
        for _ in range(node_count):
            if pos + 28 > len(raw):
                raise ValueError("truncated FEMR node data")
            node_id, x, y, z = struct.unpack_from("<iddd", raw, pos)
            pos += 28
            if node_id in nodes:
                raise ValueError(f"duplicate FEMR node id: {node_id}")
            nodes[node_id] = (x, y, z)

        elements: dict[int, FemrElement] = {}
        for _ in range(element_count):
            if pos + 6 > len(raw):
                raise ValueError("truncated FEMR element header")
            element_id, type_size = struct.unpack_from("<iH", raw, pos)
            pos += 6
            if pos + type_size + 2 > len(raw):
                raise ValueError("truncated FEMR element type")
            element_type = raw[pos:pos + type_size].decode("utf-8")
            pos += type_size
            element_node_count = struct.unpack_from("<H", raw, pos)[0]
            pos += 2
            connectivity_size = element_node_count * 4
            if pos + connectivity_size > len(raw):
                raise ValueError("truncated FEMR element connectivity")
            node_ids = struct.unpack_from(f"<{element_node_count}i", raw, pos)
            pos += connectivity_size
            if element_id in elements:
                raise ValueError(f"duplicate FEMR element id: {element_id}")
            elements[element_id] = FemrElement(element_id, element_type, node_ids)

        if pos != len(raw):
            raise ValueError("unexpected trailing bytes in FEMR mesh chunk")
        return FemrMesh(nodes, elements)


def open_results(path: str | Path) -> MeshFemrResults:
    return MeshFemrResults(path)
