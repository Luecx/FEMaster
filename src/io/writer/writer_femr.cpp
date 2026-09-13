#include "writer_femr.h"

#include "../../core/logging.h"
#include "../../model/element/element.h"
#include "../../model/model_data.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <lz4.h>
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

Bytes encode_lz4(const Bytes& input) {
    logging::error(input.size() <= static_cast<std::size_t>(LZ4_MAX_INPUT_SIZE),
                   "FemrWriter: field chunk exceeds the LZ4 block size limit");
    const int source_size = static_cast<int>(input.size());
    Bytes out(static_cast<std::size_t>(LZ4_compressBound(source_size)));
    const int stored_size = LZ4_compress_default(
        reinterpret_cast<const char*>(input.data()),
        reinterpret_cast<char*>(out.data()), source_size, static_cast<int>(out.size()));
    logging::error(stored_size > 0, "FemrWriter: LZ4 compression failed");
    out.resize(static_cast<std::size_t>(stored_size));
    return out;
}

Bytes compress(const Bytes& input, FemrCompression method) {
    if (method == FemrCompression::Lz4) return encode_lz4(input);
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

std::uint8_t field_domain_code(model::FieldDomain domain) {
    switch (domain) {
        case model::FieldDomain::UNKNOWN: return 0;
        case model::FieldDomain::NODE: return 1;
        case model::FieldDomain::ELEMENT: return 2;
        default:
            logging::error(false, "FemrWriter: unsupported FEMR v1 field domain");
            return 0;
    }
}

std::uint32_t spectral_frame(const std::string& field_name) {
    const auto separator = field_name.rfind('_');
    if (separator == std::string::npos || separator + 1 == field_name.size()) return 0;
    const std::string suffix = field_name.substr(separator + 1);
    if (!std::all_of(suffix.begin(), suffix.end(), [](unsigned char c) { return std::isdigit(c); })) return 0;
    const int one_based = std::stoi(suffix);
    return one_based > 0 ? static_cast<std::uint32_t>(one_based - 1) : 0;
}

} // namespace

FemrWriter::FemrWriter(const std::string& filename) {
    if (!filename.empty()) open(filename);
}

FemrWriter::~FemrWriter() { close(); }

FemrWriter::FemrWriter(FemrWriter&& other) noexcept
    : file_(std::move(other.file_)),
      current_loadcase_(other.current_loadcase_), current_step_type_(other.current_step_type_),
      current_frame_(other.current_frame_),
      next_field_id_(other.next_field_id_), last_frame_value_(other.last_frame_value_),
      frame_written_(other.frame_written_), frames_written_(std::move(other.frames_written_)),
      closed_(other.closed_),
      file_crc_state_(other.file_crc_state_) { other.closed_ = true; }

FemrWriter& FemrWriter::operator=(FemrWriter&& other) noexcept {
    if (this != &other) {
        close(); file_ = std::move(other.file_);
        current_loadcase_ = other.current_loadcase_; current_step_type_ = other.current_step_type_;
        current_frame_ = other.current_frame_;
        next_field_id_ = other.next_field_id_; last_frame_value_ = other.last_frame_value_;
        frame_written_ = other.frame_written_; frames_written_ = std::move(other.frames_written_);
        closed_ = other.closed_;
        file_crc_state_ = other.file_crc_state_; other.closed_ = true;
    }
    return *this;
}

void FemrWriter::open(const std::string& filename) {
    close();
    file_.open(filename, std::ios::binary | std::ios::trunc);
    logging::error(file_.is_open(), "FemrWriter: failed to open file: ", filename);
    closed_ = false; file_crc_state_ = 0xffffffffu;
    current_loadcase_ = 0; current_step_type_ = WriterStepType::Static;
    current_frame_ = 0; next_field_id_ = 1;
    frame_written_ = false; frames_written_.clear();
    last_frame_value_ = std::numeric_limits<double>::quiet_NaN();
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
    payload.push_back(static_cast<std::uint8_t>(FemrCompression::Lz4));
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
    Bytes instances;
    append_uint(instances, static_cast<std::uint32_t>(model_data.instances.size()));
    for (const auto& [name, instance] : model_data.instances) {
        logging::error(instance != nullptr, "FemrWriter: null instance: ", name);
        append_i32(instances, static_cast<std::int32_t>(instance->instance_id));
        append_string(instances, instance->instance_id == 0 ? std::string{} : name);
    }
    write_chunk("INST", instances);

    Bytes payload;
    const std::uint64_t node_count = model_data.positions
        ? static_cast<std::uint64_t>(model_data.positions->rows) : 0;
    const std::uint64_t element_count = static_cast<std::uint64_t>(model_data.elements.size());
    append_uint(payload, node_count); append_uint(payload, element_count);
    logging::error(model_data.positions != nullptr || node_count == 0,
                   "FemrWriter: mesh nodes require a positions field");
    logging::error(model_data.node_mapping.size() == static_cast<std::size_t>(node_count),
                   "FemrWriter: node mapping does not match NODE rows");
    logging::error(model_data.element_mapping.size() == static_cast<std::size_t>(element_count),
                   "FemrWriter: element mapping does not match ELEMENT rows");
    // MESH order is the normative dense field-row mapping: node/element entry
    // zero corresponds to row zero in NODE/ELEMENT FDAT arrays. Identifiers are
    // stored explicitly as well, so readers never infer them from connectivity.
    for (Index row = 0; row < static_cast<Index>(node_count); ++row) {
        const auto& [instance, local_id] = model_data.node_mapping[static_cast<std::size_t>(row)];
        logging::error(instance != nullptr, "FemrWriter: dense node has no semantic identity: ", row);
        append_i32(payload, static_cast<std::int32_t>(row));
        append_i32(payload, static_cast<std::int32_t>(instance->instance_id));
        append_i32(payload, static_cast<std::int32_t>(local_id));
        for (Index component = 0; component < 3; ++component)
            append_f64(payload, static_cast<double>((*model_data.positions)(row, component)));
    }
    for (Index row = 0; row < static_cast<Index>(element_count); ++row) {
        const auto& element = model_data.elements[static_cast<std::size_t>(row)];
        logging::error(element != nullptr, "FemrWriter: dense element row is empty: ", row);
        logging::error(element->elem_id == static_cast<ID>(row),
                       "FemrWriter: element identifier does not match its dense row");
        const auto& [instance, local_id] = model_data.element_mapping[static_cast<std::size_t>(row)];
        logging::error(instance != nullptr, "FemrWriter: dense element has no semantic identity: ", row);
        append_i32(payload, element->elem_id);
        append_i32(payload, static_cast<std::int32_t>(instance->instance_id));
        append_i32(payload, static_cast<std::int32_t>(local_id));
        append_string(payload, element->type_name());
        append_uint(payload, static_cast<std::uint16_t>(element->n_nodes()));
        for (Dim i = 0; i < element->n_nodes(); ++i) append_i32(payload, element->nodes()[i]);
    }
    write_chunk("MESH", payload);
}

void FemrWriter::add_loadcase(int id, WriterStepType step_type) {
    current_loadcase_ = id; current_step_type_ = step_type;
    current_frame_ = 0; frame_written_ = false; frames_written_.clear();
    last_frame_value_ = std::numeric_limits<double>::quiet_NaN();
    Bytes payload; append_i32(payload, id); payload.push_back(step_code(step_type));
    payload.insert(payload.end(), 3, 0);
    write_chunk("LCAS", payload);
}

void FemrWriter::ensure_frame(const std::string& field_name, Precision frame_value) {
    const double value = static_cast<double>(frame_value);
    std::uint32_t frame = current_frame_;
    if (current_step_type_ == WriterStepType::Eigenfrequency
        || current_step_type_ == WriterStepType::Buckling) {
        frame = spectral_frame(field_name);
    } else if (frame_written_ && std::isfinite(value)
               && (!std::isfinite(last_frame_value_) || value != last_frame_value_)) {
        frame = current_frame_ + 1;
    }
    if (frame_written_ && frame == current_frame_) return;
    current_frame_ = frame;
    if (frames_written_.count(frame) != 0) return;
    Bytes payload; append_i32(payload, current_loadcase_); append_uint(payload, current_frame_);
    append_f64(payload, value);
    write_chunk("FRAM", payload);
    frames_written_.insert(frame);
    last_frame_value_ = value; frame_written_ = true;
}

void FemrWriter::write_field(const model::Field& field, const std::string& field_name,
                             const model::ModelData*, Precision frame_value) {
    logging::error(file_.is_open(), "FemrWriter: file is not open");
    const bool supported_domain = field.domain == model::FieldDomain::UNKNOWN
                               || field.domain == model::FieldDomain::NODE
                               || field.domain == model::FieldDomain::ELEMENT;
    if (!supported_domain) {
        logging::warning(false,
                         "FemrWriter v1 skips unsupported field domain for field: ",
                         field_name);
        return;
    }
    ensure_frame(field_name, frame_value);
    const std::uint64_t field_id = next_field_id_++;
    Bytes metadata;
    append_uint(metadata, field_id); append_i32(metadata, current_loadcase_);
    append_uint(metadata, current_frame_); append_string(metadata, field_name);
    metadata.push_back(field_domain_code(field.domain));
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
    write_chunk("FDAT", data, FemrCompression::Lz4);
}

} // namespace writer
} // namespace io
} // namespace fem
