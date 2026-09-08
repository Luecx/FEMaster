# FEMR binary result format (version 1.0)

FEMR is FEMaster's seekable binary result container. Multi-byte integers and
floating-point values are little-endian. Unknown chunks must be skipped using
their stored size, which permits compatible format extensions.


Binary results can be emitted alongside the existing `.res` and `.frd` files:

```text
FEMaster model.inp --result-format binary --result-compression zstd
```

The resulting `.femr` container is chunked and supports lazy Python access via
`femaster_api.open_results("model.femr")`. See
[`documentation/FEMR_FORMAT.md`](documentation/FEMR_FORMAT.md) for the versioned
on-disk specification.

## Chunk envelope

Every chunk starts with the same 32-byte envelope:

| Offset | Type | Meaning |
|---:|---|---|
| 0 | `char[4]` | chunk code |
| 4 | `uint32` | compression: 0 none, 1 LZ4 block, 2 Zstandard frame |
| 8 | `uint64` | stored payload bytes |
| 16 | `uint64` | uncompressed payload bytes |
| 24 | `uint32` | CRC32 of the uncompressed payload |
| 28 | `uint32` | reserved, zero in version 1 |

Only `FDAT` is compressed in version 1. Each field is compressed separately,
so reading one field never requires decompressing another field.

## Chunks

- `HEAD` — must be first. Payload: magic `FEMR`, major/minor `uint16`, byte
  order (1 = little), scalar width (4 or 8), default compression, reserved,
  and a reserved `uint64`.
- `MESH` — node and element topology. Starts with node/element counts as
  `uint64`. A node is `int32 id` plus three `float64` coordinates. An element
  is `int32 id`, a length-prefixed UTF-8 type name, `uint16` node count, and
  `int32` node identifiers.
- `LCAS` — `int32 id`, `uint8` step type, three reserved bytes. Step types are
  static (0), dynamic (1), eigenfrequency (2), and buckling (3).
- `FRAM` — `int32 loadcase id`, `uint32 frame id`, `float64 frame value`.
- `FMET` — `uint64 field id`, loadcase/frame ids, a length-prefixed UTF-8 name,
  domain and scalar-width bytes, two reserved bytes, then row and component
  counts as `uint64`.
- `FDAT` — `uint64 field id` followed by the dense row-major scalar array.
- `CSUM` — must be last. Its `uint32` payload is the CRC32 of every preceding
  chunk envelope and stored payload byte.

Strings use a `uint16` byte length followed by UTF-8 bytes. Field domains are
unknown (0), node (1), element (2), element-nodal (3), element-IP (4), and
element-material-point (5).

## Lazy access

A reader scans only envelopes plus `HEAD`, `LCAS`, `FRAM`, and `FMET`. It skips
`MESH` and every `FDAT` by seeking over `stored payload bytes`, retaining each
field-data offset. The field's `FDAT` payload is read, decompressed, and CRC
checked only when requested.

Python usage:

```python
import femaster_api as femaster

with femaster.open_results("model.femr") as results:
    stress = results.loadcase(1).frame(20).field("STRESS")
```
