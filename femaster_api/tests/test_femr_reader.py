from __future__ import annotations

import struct
import zlib

import femaster_api as femaster


CHUNK = struct.Struct("<4sIQQII")


def _chunk(kind: bytes, payload: bytes) -> bytes:
    return CHUNK.pack(kind, 0, len(payload), len(payload), zlib.crc32(payload), 0) + payload


def test_open_results_indexes_then_lazily_reads_one_field(tmp_path) -> None:
    head = b"FEMR" + struct.pack("<HHBBBBQ", 1, 0, 1, 8, 0, 0, 0)
    loadcase = struct.pack("<iB3x", 1, 0)
    frame = struct.pack("<iId", 1, 20, 2.0)
    name = b"STRESS"
    metadata = struct.pack("<QiIH", 7, 1, 20, len(name)) + name
    metadata += struct.pack("<BBHQQ", 1, 8, 0, 2, 3)
    values = struct.pack("<Q6d", 7, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    path = tmp_path / "model.femr"
    path.write_bytes(
        _chunk(b"HEAD", head)
        + _chunk(b"LCAS", loadcase)
        + _chunk(b"FRAM", frame)
        + _chunk(b"FMET", metadata)
        + _chunk(b"FDAT", values)
        + _chunk(b"CSUM", struct.pack("<I", 0))
    )

    with femaster.open_results(path) as results:
        lazy_frame = results.loadcase(1).frame(20)
        assert lazy_frame.field_names == ("STRESS",)
        assert lazy_frame._cache == {}
        stress = lazy_frame.field("STRESS")
        assert stress.row(1) == (4.0, 5.0, 6.0)
        assert lazy_frame.field("STRESS") is stress
