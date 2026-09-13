#pragma once

#include "../../core/core.h"
#include "../../data/field.h"
#include "writer_step_type.h"

#include <cstdint>
#include <fstream>
#include <limits>
#include <set>
#include <string>
#include <vector>

namespace fem {
namespace model { struct ModelData; }
namespace io {
namespace writer {

enum class FemrCompression : std::uint8_t { None = 0, Lz4 = 1 };

/** Versioned, chunked FEMaster binary result writer.
 *
 * Field metadata and field bytes are separate chunks.  A reader can scan chunk
 * headers and FMET payloads once, then seek directly to an FDAT chunk without
 * reading unrelated arrays.
 */
class FemrWriter {
public:
    explicit FemrWriter(const std::string& filename = "");
    ~FemrWriter();

    FemrWriter(FemrWriter&& other) noexcept;
    FemrWriter& operator=(FemrWriter&& other) noexcept;
    FemrWriter(const FemrWriter&) = delete;
    FemrWriter& operator=(const FemrWriter&) = delete;

    void open(const std::string& filename);
    void close();
    void write_model_data(const model::ModelData& model_data);
    void add_loadcase(int id, WriterStepType step_type = WriterStepType::Static);
    void write_field(const model::Field& field,
                     const std::string& field_name,
                     const model::ModelData* model_data = nullptr,
                     Precision frame_value = std::numeric_limits<Precision>::quiet_NaN());

private:
    void write_header();
    void write_chunk(const char type[4], const std::vector<std::uint8_t>& payload,
                     FemrCompression compression = FemrCompression::None,
                     bool include_in_file_checksum = true);
    void ensure_frame(const std::string& field_name, Precision frame_value);

    std::ofstream file_;
    int current_loadcase_{0};
    WriterStepType current_step_type_{WriterStepType::Static};
    std::uint32_t current_frame_{0};
    std::uint64_t next_field_id_{1};
    double last_frame_value_{std::numeric_limits<double>::quiet_NaN()};
    bool frame_written_{false};
    std::set<std::uint32_t> frames_written_;
    bool closed_{true};
    std::uint32_t file_crc_state_{0xffffffffu};
};

} // namespace writer
} // namespace io
} // namespace fem
