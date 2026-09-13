"""Lazy mesh extension for the FEMR result reader."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import struct

from .femr_reader import FemrResults, _CHUNK, _DataRef, _decompress, _verify


@dataclass(frozen=True, slots=True)
class FemrElement:
    id: int
    instance_id: int
    local_id: int
    type: str
    node_ids: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class FemrMesh:
    nodes: dict[int, tuple[float, float, float]]
    elements: dict[int, FemrElement]
    node_ids: tuple[int, ...]
    element_ids: tuple[int, ...]
    node_semantic_ids: dict[int, tuple[int, int]]
    instance_names: dict[int, str]

    def node(self, id: int) -> tuple[float, float, float]:
        return self.nodes[id]

    def element(self, id: int) -> FemrElement:
        return self.elements[id]

    def node_at(self, row: int) -> tuple[float, float, float]:
        """Return the node represented by a NODE field row."""
        return self.nodes[self.node_ids[row]]

    def element_at(self, row: int) -> FemrElement:
        """Return the element represented by an ELEMENT field row."""
        return self.elements[self.element_ids[row]]

    def node_label(self, id: int) -> str:
        instance_id, local_id = self.node_semantic_ids[id]
        name = self.instance_names.get(instance_id, "")
        return str(local_id) if not name else f"{name}.{local_id}"

    def element_label(self, id: int) -> str:
        element = self.elements[id]
        name = self.instance_names.get(element.instance_id, "")
        return str(element.local_id) if not name else f"{name}.{element.local_id}"


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
        node_ids: list[int] = []
        node_semantic_ids: dict[int, tuple[int, int]] = {}
        for _ in range(node_count):
            if pos + 36 > len(raw):
                raise ValueError("truncated FEMR node data")
            node_id, instance_id, local_id, x, y, z = struct.unpack_from("<iiiddd", raw, pos)
            pos += 36
            if node_id in nodes:
                raise ValueError(f"duplicate FEMR node id: {node_id}")
            if instance_id not in self._instances:
                raise ValueError(f"unknown FEMR node instance id: {instance_id}")
            nodes[node_id] = (x, y, z)
            node_ids.append(node_id)
            node_semantic_ids[node_id] = (instance_id, local_id)

        elements: dict[int, FemrElement] = {}
        element_ids: list[int] = []
        for _ in range(element_count):
            if pos + 14 > len(raw):
                raise ValueError("truncated FEMR element header")
            element_id, instance_id, local_id, type_size = struct.unpack_from("<iiiH", raw, pos)
            pos += 14
            if instance_id not in self._instances:
                raise ValueError(f"unknown FEMR element instance id: {instance_id}")
            if pos + type_size + 2 > len(raw):
                raise ValueError("truncated FEMR element type")
            element_type = raw[pos:pos + type_size].decode("utf-8")
            pos += type_size
            element_node_count = struct.unpack_from("<H", raw, pos)[0]
            pos += 2
            connectivity_size = element_node_count * 4
            if pos + connectivity_size > len(raw):
                raise ValueError("truncated FEMR element connectivity")
            connectivity = struct.unpack_from(f"<{element_node_count}i", raw, pos)
            pos += connectivity_size
            if element_id in elements:
                raise ValueError(f"duplicate FEMR element id: {element_id}")
            elements[element_id] = FemrElement(
                element_id, instance_id, local_id, element_type, connectivity
            )
            element_ids.append(element_id)

        if pos != len(raw):
            raise ValueError("unexpected trailing bytes in FEMR mesh chunk")
        return FemrMesh(
            nodes, elements, tuple(node_ids), tuple(element_ids),
            node_semantic_ids, self._instances
        )


def open_results(path: str | Path) -> MeshFemrResults:
    return MeshFemrResults(path)
