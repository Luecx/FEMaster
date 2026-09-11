"""Lazy reader for chunked FEMaster ``.femr`` result files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import struct
import zlib

from femaster_api.model.fields import Field, FieldDomain

_CHUNK = struct.Struct("<4sIQQII")
_DOMAIN = {
    0: FieldDomain.UNKNOWN,
    1: FieldDomain.NODE,
    2: FieldDomain.ELEMENT,
}


@dataclass(frozen=True, slots=True)
class _DataRef:
    offset: int
    stored_size: int
    raw_size: int
    compression: int
    checksum: int


@dataclass(frozen=True, slots=True)
class _FieldMeta:
    field_id: int
    loadcase_id: int
    frame_id: int
    name: str
    domain: int
    dtype_size: int
    rows: int
    components: int


class FemrFrame:
    """A frame whose field arrays remain on disk until requested."""

    def __init__(self, owner: "FemrResults", id: int, value: float) -> None:
        self._owner = owner
        self.id = id
        self.value = value
        self._fields: dict[str, _FieldMeta] = {}
        self._cache: dict[str, Field] = {}

    @property
    def field_names(self) -> tuple[str, ...]:
        return tuple(self._fields)

    def field(self, name: str) -> Field:
        """Load and cache one field without touching other field-data chunks."""
        if name not in self._cache:
            self._cache[name] = self._owner._load_field(self._fields[name])
        return self._cache[name]


class FemrLoadCase:
    def __init__(self, id: int, step_type: int) -> None:
        self.id = id
        self.step_type = step_type
        self._frames: dict[int, FemrFrame] = {}

    @property
    def frame_ids(self) -> tuple[int, ...]:
        return tuple(self._frames)

    def frame(self, id: int = 0) -> FemrFrame:
        return self._frames[id]


class FemrResults:
    """Seekable FEMR index. Use as a context manager for deterministic close."""

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self._file = self.path.open("rb")
        self._loadcases: dict[int, FemrLoadCase] = {}
        self._data: dict[int, _DataRef] = {}
        try:
            self._scan()
        except Exception:
            self._file.close()
            raise

    @property
    def loadcase_ids(self) -> tuple[int, ...]:
        return tuple(self._loadcases)

    def loadcase(self, id: int = 1) -> FemrLoadCase:
        return self._loadcases[id]

    def close(self) -> None:
        self._file.close()

    def __enter__(self) -> "FemrResults":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def _scan(self) -> None:
        current_loadcase: FemrLoadCase | None = None
        current_frame: FemrFrame | None = None
        metadata: dict[int, _FieldMeta] = {}
        pending_field_id: int | None = None
        saw_header = False

        while True:
            raw_header = self._file.read(_CHUNK.size)
            if not raw_header:
                break
            if len(raw_header) != _CHUNK.size:
                raise ValueError("truncated FEMR chunk header")
            kind, compression, stored_size, raw_size, checksum, _ = _CHUNK.unpack(raw_header)
            payload_offset = self._file.tell()

            if kind == b"FDAT":
                if pending_field_id is None:
                    raise ValueError("FEMR field-data chunk has no preceding metadata")
                self._data[pending_field_id] = _DataRef(
                    payload_offset, stored_size, raw_size, compression, checksum
                )
                pending_field_id = None
                self._file.seek(stored_size, 1)
                continue

            if kind == b"MESH":
                # Mesh access can be added independently; opening results must
                # not allocate or checksum a potentially very large mesh.
                self._file.seek(stored_size, 1)
                continue

            stored = self._file.read(stored_size)
            raw = _decompress(stored, compression, raw_size)
            _verify(raw, checksum, kind)
            if kind == b"HEAD":
                if len(raw) < 12 or raw[:4] != b"FEMR":
                    raise ValueError("not a FEMR file")
                major, minor = struct.unpack_from("<HH", raw, 4)
                if major != 1:
                    raise ValueError(f"unsupported FEMR version {major}.{minor}")
                if raw[8] != 1:
                    raise ValueError("unsupported FEMR byte order")
                saw_header = True
            elif kind == b"LCAS":
                loadcase_id, step_type = struct.unpack_from("<iB", raw)
                current_loadcase = FemrLoadCase(loadcase_id, step_type)
                self._loadcases[loadcase_id] = current_loadcase
                current_frame = None
            elif kind == b"FRAM":
                loadcase_id, frame_id, value = struct.unpack_from("<iId", raw)
                current_loadcase = self._loadcases.setdefault(loadcase_id, FemrLoadCase(loadcase_id, 0))
                current_frame = FemrFrame(self, frame_id, value)
                current_loadcase._frames[frame_id] = current_frame
            elif kind == b"FMET":
                meta = _parse_metadata(raw)
                metadata[meta.field_id] = meta
                pending_field_id = meta.field_id
                loadcase = self._loadcases.setdefault(meta.loadcase_id, FemrLoadCase(meta.loadcase_id, 0))
                frame = loadcase._frames.setdefault(meta.frame_id, FemrFrame(self, meta.frame_id, float("nan")))
                frame._fields[meta.name] = meta

        if not saw_header:
            raise ValueError("missing FEMR header chunk")
        missing = metadata.keys() - self._data.keys()
        if missing:
            raise ValueError(f"missing FEMR field-data chunks: {sorted(missing)}")

    def _load_field(self, meta: _FieldMeta) -> Field:
        if meta.domain not in _DOMAIN:
            raise ValueError(f"unsupported FEMR v1 field domain: {meta.domain}")
        ref = self._data[meta.field_id]
        self._file.seek(ref.offset)
        stored = self._file.read(ref.stored_size)
        raw = _decompress(stored, ref.compression, ref.raw_size)
        _verify(raw, ref.checksum, b"FDAT")
        if struct.unpack_from("<Q", raw)[0] != meta.field_id:
            raise ValueError("FEMR field metadata/data id mismatch")
        fmt = "f" if meta.dtype_size == 4 else "d" if meta.dtype_size == 8 else None
        if fmt is None:
            raise ValueError(f"unsupported FEMR scalar width: {meta.dtype_size}")
        count = meta.rows * meta.components
        expected = 8 + count * meta.dtype_size
        if len(raw) != expected:
            raise ValueError(f"FEMR field {meta.name!r} has invalid data size")
        flat = struct.unpack_from(f"<{count}{fmt}", raw, 8)
        values = {row: tuple(flat[row * meta.components:(row + 1) * meta.components]) for row in range(meta.rows)}
        return Field(meta.name, _DOMAIN.get(meta.domain, FieldDomain.UNKNOWN), meta.components, values, "NONE")


def open_results(path: str | Path) -> FemrResults:
    """Open a ``.femr`` file and build its lightweight lazy index."""
    return FemrResults(path)


def _parse_metadata(raw: bytes) -> _FieldMeta:
    if len(raw) < 28:
        raise ValueError("invalid FEMR field metadata")
    field_id, loadcase_id, frame_id, name_size = struct.unpack_from("<QiIH", raw)
    pos = 18
    end = pos + name_size
    if end + 20 > len(raw):
        raise ValueError("truncated FEMR field metadata")
    name = raw[pos:end].decode("utf-8")
    domain, dtype_size = struct.unpack_from("<BB", raw, end)
    rows, components = struct.unpack_from("<QQ", raw, end + 4)
    return _FieldMeta(field_id, loadcase_id, frame_id, name, domain, dtype_size, rows, components)


def _verify(raw: bytes, expected: int, kind: bytes) -> None:
    if zlib.crc32(raw) & 0xFFFFFFFF != expected:
        raise ValueError(f"FEMR checksum mismatch in {kind.decode('ascii', 'replace')} chunk")


def _decompress(stored: bytes, method: int, raw_size: int) -> bytes:
    if method == 0:
        raw = stored
    elif method == 1:
        raw = _decode_lz4(stored)
    else:
        raise ValueError(f"unsupported FEMR compression id: {method}")
    if len(raw) != raw_size:
        raise ValueError("FEMR decompressed size mismatch")
    return raw


def _decode_lz4(data: bytes) -> bytes:
    out = bytearray(); pos = 0
    while pos < len(data):
        token = data[pos]; pos += 1
        literal = token >> 4
        if literal == 15:
            while True:
                value = data[pos]; pos += 1; literal += value
                if value != 255: break
        out += data[pos:pos + literal]; pos += literal
        if pos == len(data): break
        offset = data[pos] | (data[pos + 1] << 8); pos += 2
        if not offset or offset > len(out): raise ValueError("invalid LZ4 offset")
        match = token & 15
        if match == 15:
            while True:
                value = data[pos]; pos += 1; match += value
                if value != 255: break
        match += 4
        for _ in range(match): out.append(out[-offset])
    return bytes(out)
