#include "writer_femr.h"

#include "../../core/logging.h"
#include "../../model/element/element.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <set>
#include <utility>

namespace fem {
namespace io {
namespace writer {
namespace {

using Bytes = std::vector<std::uint8_t>;

std::uint32_t crc32_update(std::uint32_t crc, const std::uint8_t* data, std::size_t size) {
    for (std::size_t i = 0; i < size; ++i) {
        crc ^= data[i];
        for (int bit = 0; bit < 8; ++bit)
            crc = (crc >> 1) ^ (0xedb88320u & (0u - (crc & 1u)));
    }
    return crc;
}

std::uint32_t crc32(const Bytes& data) {
    return crc32_update(0xffffffffu, data.data(), data.size()) ^ 0xffffffffu;
}

template<class UInt>
void append_uint(Bytes& out, UInt value) {
    for (std::size_t i = 0; i < sizeof(UInt); ++i)
        out.push_back(static_cast<std::uint8_t>((value >> (8 * i)) & 0xffu));
}

void append_i32(Bytes& out, std::int32_t value) {
    append_uint(out, static_cast<std::uint32_t>(value));
}

void append_f64(Bytes& out, double value) {
    std::uint64_t bits{};
    std::memcpy(&bits, &value, sizeof(bits));
    append_uint(out, bits);
}

void append_precision(Bytes& out, Precision value) {
    if (sizeof(Precision) == sizeof(double)) {
        append_f64(out, static_cast<double>(value));
    } else {
        const float as_float = static_cast<float>(value);
        std::uint32_t bits{};
        std::memcpy(&bits, &as_float, sizeof(bits));
        append_uint(out, bits);
    }
}

void append_string(Bytes& out, const std::string& value) {
    logging::error(value.size() <= 0xffffu, "FemrWriter: string is too long");
    append_uint(out, static_cast<std::uint16_t>(value.size()));
    out.insert(out.end(), value.begin(), value.end());
}

void append_lz4_length(Bytes& out, std::size_t length) {
    while (length >= 255) { out.push_back(255); length -= 255; }
    out.push_back(static_cast<std::uint8_t>(length));
}

std::uint32_t read_u32(const std::uint8_t* data) {
    return static_cast<std::uint32_t>(data[0])
         | (static_cast<std::uint32_t>(data[1]) << 8)
         | (static_cast<std::uint32_t>(data[2]) << 16)
         | (static_cast<std::uint32_t>(data[3]) << 24);
}

// Dependency-free LZ4 block compressor. Fields are independent blocks, which
// both bounds memory use and is essential for field-level lazy loading.
Bytes encode_lz4(const Bytes& input) {
    Bytes out;
    constexpr std::size_t hash_size = 1u << 16;
    std::vector<std::int64_t> table(hash_size, -1);
    std::size_t anchor = 0, cursor = 0;

    while (cursor + 4 <= input.size()) {
        const std::uint32_t sequence = read_u32(input.data() + cursor);
        const std::size_t hash = (sequence * 2654435761u) >> 16;
        const std::int64_t candidate = table[hash];
        table[hash] = static_cast<std::int64_t>(cursor);
        if (candidate < 0 || cursor - static_cast<std::size_t>(candidate) > 65535 ||
            read_u32(input.data() + candidate) != sequence) {
            ++cursor;
            continue;
        }

        const std::size_t literal = cursor - anchor;
        std::size_t match = 4;
        while (cursor + match < input.size() &&
               input[static_cast<std::size_t>(candidate) + match] == input[cursor + match]) ++match;
        const std::size_t encoded_match = match - 4;
        out.push_back(static_cast<std::uint8_t>((std::min<std::size_t>(literal, 15) << 4)
                                               | std::min<std::size_t>(encoded_match, 15)));
        if (literal >= 15) append_lz4_length(out, literal - 15);
        out.insert(out.end(), input.begin() + static_cast<std::ptrdiff_t>(anchor),
                   input.begin() + static_cast<std::ptrdiff_t>(cursor));
        const std::size_t offset = cursor - static_cast<std::size_t>(candidate);
        append_uint(out, static_cast<std::uint16_t>(offset));
        if (encoded_match >= 15) append_lz4_length(out, encoded_match - 15);
        cursor += match; anchor = cursor;
    }

    const std::size_t literal = input.size() - anchor;
    out.push_back(static_cast<std::uint8_t>(std::min<std::size_t>(literal, 15) << 4));
    if (literal >= 15) append_lz4_length(out, literal - 15);
    out.insert(out.end(), input.begin() + static_cast<std::ptrdiff_t>(anchor), input.end());
    return out;
}

// Minimal standards-compliant Zstandard encoder using raw and RLE blocks. It
// has no linked dependency and compresses repeated byte runs (common in sparse
// and zero-valued result arrays); arbitrary data remains an interoperable raw
// block within the Zstandard frame.
Bytes encode_zstd(const Bytes& input) {
    Bytes out{0x28, 0xb5, 0x2f, 0xfd};
    const std::uint64_t size = input.size();
    if (size < 256) {
        out.push_back(0x20); append_uint(out, static_cast<std::uint8_t>(size));
    } else if (size < 65792) {
        out.push_back(0x60); append_uint(out, static_cast<std::uint16_t>(size - 256));
    } else if (size <= 0xffffffffu) {
        out.push_back(0xa0); append_uint(out, static_cast<std::uint32_t>(size));
    } else {
        out.push_back(0xe0); append_uint(out, size);
    }

    struct Block { std::size_t offset; std::size_t size; bool rle; };
    constexpr std::size_t max_block = 131071, min_rle = 8;
    std::vector<Block> blocks;
    std::size_t cursor = 0;
    while (cursor < input.size()) {
        std::size_t run = 1;
        while (cursor + run < input.size() && input[cursor + run] == input[cursor] && run < max_block) ++run;
        if (run >= min_rle) {
            blocks.push_back({cursor, run, true}); cursor += run; continue;
        }
        const std::size_t raw_begin = cursor++;
        while (cursor < input.size() && cursor - raw_begin < max_block) {
            run = 1;
            while (cursor + run < input.size() && input[cursor + run] == input[cursor] && run < max_block) ++run;
            if (run >= min_rle) break;
            ++cursor;
        }
        blocks.push_back({raw_begin, cursor - raw_begin, false});
    }
    if (blocks.empty()) blocks.push_back({0, 0, false});
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        const Block block = blocks[i];
        const bool last = i + 1 == blocks.size();
        const std::uint32_t block_header = (static_cast<std::uint32_t>(block.size) << 3)
                                         | (block.rle ? 2u : 0u) | (last ? 1u : 0u);
        out.push_back(static_cast<std::uint8_t>(block_header));
        out.push_back(static_cast<std::uint8_t>(block_header >> 8));
        out.push_back(static_cast<std::uint8_t>(block_header >> 16));
        if (block.rle) out.push_back(input[block.offset]);
        else out.insert(out.end(), input.begin() + static_cast<std::ptrdiff_t>(block.offset),
                        input.begin() + static_cast<std::ptrdiff_t>(block.offset + block.size));
    }
    return out;
}

Bytes compress(const Bytes& input, FemrCompression method) {
    if (method == FemrCompression::Lz4) return encode_lz4(input);
    if (method == FemrCompression::Zstd) return encode_zstd(input);
    return input;
}

std::uint8_t step_code(WriterStepType type) {
    switch (type) {
        case WriterStepType::Static: return 0;
        case WriterStepType::Dynamic: return 1;
        case WriterStepType::Eigenfrequency: return 2;
        case WriterStepType::Buckling: return 3;
    }
    return 0;
}

} // namespace

FemrCompression femr_compression_from_string(const std::string& value) {
    if (value == "none") return FemrCompression::None;
    if (value == "lz4") return FemrCompression::Lz4;
    if (value == "zstd") return FemrCompression::Zstd;
    logging::error(false, "FemrWriter: unknown compression: ", value);
    return FemrCompression::None;
}

FemrWriter::FemrWriter(const std::string& filename, FemrCompression compression)
    : compression_(compression) {
    if (!filename.empty()) open(filename);
}

FemrWriter::~FemrWriter() { close(); }

FemrWriter::FemrWriter(FemrWriter&& other) noexcept
    : file_(std::move(other.file_)), compression_(other.compression_),
      current_loadcase_(other.current_loadcase_), current_frame_(other.current_frame_),
      next_field_id_(other.next_field_id_), last_frame_value_(other.last_frame_value_),
      frame_written_(other.frame_written_), closed_(other.closed_),
      file_crc_state_(other.file_crc_state_) { other.closed_ = true; }

FemrWriter& FemrWriter::operator=(FemrWriter&& other) noexcept {
    if (this != &other) {
        close(); file_ = std::move(other.file_); compression_ = other.compression_;
        current_loadcase_ = other.current_loadcase_; current_frame_ = other.current_frame_;
        next_field_id_ = other.next_field_id_; last_frame_value_ = other.last_frame_value_;
        frame_written_ = other.frame_written_; closed_ = other.closed_;
        file_crc_state_ = other.file_crc_state_; other.closed_ = true;
    }
    return *this;
}

void FemrWriter::open(const std::string& filename) {
    close();
    file_.open(filename, std::ios::binary | std::ios::trunc);
    logging::error(file_.is_open(), "FemrWriter: failed to open file: ", filename);
    closed_ = false; file_crc_state_ = 0xffffffffu;
    current_loadcase_ = 0; current_frame_ = 0; next_field_id_ = 1;
    frame_written_ = false; last_frame_value_ = std::numeric_limits<double>::quiet_NaN();
    write_header();
}

void FemrWriter::close() {
    if (!file_.is_open() || closed_) return;
    Bytes checksum; append_uint(checksum, file_crc_state_ ^ 0xffffffffu);
    write_chunk("CSUM", checksum, FemrCompression::None, false);
    file_.close(); closed_ = true;
}

void FemrWriter::write_header() {
    Bytes payload{'F','E','M','R'};
    append_uint(payload, std::uint16_t{1});
    append_uint(payload, std::uint16_t{0});
    payload.push_back(1); // little endian
    payload.push_back(static_cast<std::uint8_t>(sizeof(Precision)));
    payload.push_back(static_cast<std::uint8_t>(compression_));
    payload.push_back(0);
    append_uint(payload, std::uint64_t{0});
    write_chunk("HEAD", payload);
}

void FemrWriter::write_chunk(const char type[4], const Bytes& payload,
                             FemrCompression compression, bool include_in_file_checksum) {
    const Bytes stored = compress(payload, compression);
    Bytes header;
    header.insert(header.end(), type, type + 4);
    append_uint(header, static_cast<std::uint32_t>(compression));
    append_uint(header, static_cast<std::uint64_t>(stored.size()));
    append_uint(header, static_cast<std::uint64_t>(payload.size()));
    append_uint(header, crc32(payload));
    append_uint(header, std::uint32_t{0});
    file_.write(reinterpret_cast<const char*>(header.data()), static_cast<std::streamsize>(header.size()));
    file_.write(reinterpret_cast<const char*>(stored.data()), static_cast<std::streamsize>(stored.size()));
    logging::error(file_.good(), "FemrWriter: failed while writing chunk");
    if (include_in_file_checksum) {
        file_crc_state_ = crc32_update(file_crc_state_, header.data(), header.size());
        file_crc_state_ = crc32_update(file_crc_state_, stored.data(), stored.size());
    }
}

void FemrWriter::write_model_data(const model::ModelData& model_data) {
    Bytes payload;
    std::set<ID> node_ids;
    for (const auto& element : model_data.elements) {
        if (!element) continue;
        for (ID node_id : *element) node_ids.insert(node_id);
    }
    const std::uint64_t node_count = node_ids.size();
    const std::uint64_t element_count = static_cast<std::uint64_t>(std::count_if(
        model_data.elements.begin(), model_data.elements.end(), [](const auto& item) { return item != nullptr; }));
    append_uint(payload, node_count); append_uint(payload, element_count);
    logging::error(model_data.positions != nullptr || node_count == 0,
                   "FemrWriter: mesh nodes require a positions field");
    for (ID id : node_ids) {
        append_i32(payload, id);
        for (Index component = 0; component < 3; ++component)
            append_f64(payload, static_cast<double>((*model_data.positions)(static_cast<Index>(id), component)));
    }
    for (const auto& element : model_data.elements) {
        if (!element) continue;
        append_i32(payload, element->elem_id);
        append_string(payload, element->type_name());
        append_uint(payload, static_cast<std::uint16_t>(element->n_nodes()));
        for (Dim i = 0; i < element->n_nodes(); ++i) append_i32(payload, element->nodes()[i]);
    }
    write_chunk("MESH", payload);
}

void FemrWriter::add_loadcase(int id, WriterStepType step_type) {
    current_loadcase_ = id; current_frame_ = 0; frame_written_ = false;
    last_frame_value_ = std::numeric_limits<double>::quiet_NaN();
    Bytes payload; append_i32(payload, id); payload.push_back(step_code(step_type));
    payload.insert(payload.end(), 3, 0);
    write_chunk("LCAS", payload);
}

void FemrWriter::ensure_frame(Precision frame_value) {
    const double value = static_cast<double>(frame_value);
    if (frame_written_ && std::isfinite(value) &&
        (!std::isfinite(last_frame_value_) || value != last_frame_value_)) {
        ++current_frame_; frame_written_ = false;
    }
    if (frame_written_) return;
    Bytes payload; append_i32(payload, current_loadcase_); append_uint(payload, current_frame_);
    append_f64(payload, value);
    write_chunk("FRAM", payload);
    last_frame_value_ = value; frame_written_ = true;
}

void FemrWriter::write_field(const model::Field& field, const std::string& field_name,
                             const model::ModelData*, Precision frame_value) {
    logging::error(file_.is_open(), "FemrWriter: file is not open");
    ensure_frame(frame_value);
    const std::uint64_t field_id = next_field_id_++;
    Bytes metadata;
    append_uint(metadata, field_id); append_i32(metadata, current_loadcase_);
    append_uint(metadata, current_frame_); append_string(metadata, field_name);
    metadata.push_back(static_cast<std::uint8_t>(field.domain));
    metadata.push_back(static_cast<std::uint8_t>(sizeof(Precision)));
    append_uint(metadata, std::uint16_t{0});
    append_uint(metadata, static_cast<std::uint64_t>(field.rows));
    append_uint(metadata, static_cast<std::uint64_t>(field.components));
    write_chunk("FMET", metadata);

    Bytes values;
    values.reserve(field.values.size() * sizeof(Precision));
    for (Index row = 0; row < field.rows; ++row)
        for (Index component = 0; component < field.components; ++component)
            append_precision(values, field(row, component));
    Bytes data; append_uint(data, field_id); data.insert(data.end(), values.begin(), values.end());
    write_chunk("FDAT", data, compression_);
}

} // namespace writer
} // namespace io
} // namespace fem
