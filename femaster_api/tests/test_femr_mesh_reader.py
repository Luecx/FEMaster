from __future__ import annotations

import struct
import zlib

import femaster_api as femaster


CHUNK = struct.Struct("<4sIQQII")


def _chunk(kind: bytes, payload: bytes) -> bytes:
    return CHUNK.pack(kind, 0, len(payload), len(payload), zlib.crc32(payload), 0) + payload


def test_mesh_is_loaded_once_and_only_on_access(tmp_path) -> None:
    head = b"FEMR" + struct.pack("<HHBBBBQ", 1, 0, 1, 8, 0, 0, 0)
    instances = struct.pack("<I", 2)
    instances += struct.pack("<iH", 0, 0)
    instances += struct.pack("<iH", 1, 4) + b"bolt"
    mesh = struct.pack("<QQ", 2, 1)
    mesh += struct.pack("<iiiddd", 10, 0, 10, 0.0, 1.0, 2.0)
    mesh += struct.pack("<iiiddd", 20, 1, 5, 3.0, 4.0, 5.0)
    mesh += struct.pack("<iiiH", 7, 1, 3, 2) + b"T3" + struct.pack("<H2i", 2, 10, 20)
    path = tmp_path / "mesh.femr"
    path.write_bytes(_chunk(b"HEAD", head) + _chunk(b"INST", instances) + _chunk(b"MESH", mesh))

    with femaster.open_results(path) as results:
        assert results._mesh is None
        loaded = results.mesh
        assert loaded.node(10) == (0.0, 1.0, 2.0)
        assert loaded.element(7).type == "T3"
        assert loaded.element(7).node_ids == (10, 20)
        assert loaded.node_ids == (10, 20)
        assert loaded.element_ids == (7,)
        assert loaded.node_at(1) == (3.0, 4.0, 5.0)
        assert loaded.element_at(0).id == 7
        assert loaded.node_label(10) == "10"
        assert loaded.node_label(20) == "bolt.5"
        assert loaded.element_label(7) == "bolt.3"
        assert results.instances == {0: "", 1: "bolt"}
        assert results.mesh is loaded
