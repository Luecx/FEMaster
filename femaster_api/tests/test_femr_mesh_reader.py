from __future__ import annotations

import struct
import zlib

import femaster_api as femaster


CHUNK = struct.Struct("<4sIQQII")


def _chunk(kind: bytes, payload: bytes) -> bytes:
    return CHUNK.pack(kind, 0, len(payload), len(payload), zlib.crc32(payload), 0) + payload


def test_mesh_is_loaded_once_and_only_on_access(tmp_path) -> None:
    head = b"FEMR" + struct.pack("<HHBBBBQ", 1, 0, 1, 8, 0, 0, 0)
    mesh = struct.pack("<QQ", 2, 1)
    mesh += struct.pack("<iddd", 10, 0.0, 1.0, 2.0)
    mesh += struct.pack("<iddd", 20, 3.0, 4.0, 5.0)
    mesh += struct.pack("<iH", 7, 2) + b"T3" + struct.pack("<H2i", 2, 10, 20)
    path = tmp_path / "mesh.femr"
    path.write_bytes(_chunk(b"HEAD", head) + _chunk(b"MESH", mesh))

    with femaster.open_results(path) as results:
        assert results._mesh is None
        loaded = results.mesh
        assert loaded.node(10) == (0.0, 1.0, 2.0)
        assert loaded.element(7).type == "T3"
        assert loaded.element(7).node_ids == (10, 20)
        assert results.mesh is loaded
