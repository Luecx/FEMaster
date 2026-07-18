#pragma once

#include "../../core/core.h"
#include "../../core/logging.h"
#include "../../data/field.h"
#include "writer_step_type.h"

#include <fstream>
#include <limits>
#include <string>

namespace fem {
namespace model {
struct ModelData;
}

namespace io {
namespace writer {

/**
 * @brief FEMaster .res writer.
 *
 * Writes all model::Field instances to the native FEMaster result format.
 * The field domain is taken from field.domain. ELEMENT_NODAL, ELEMENT_IP and ELEMENT_MP
 * fields require ModelData so that element/local-location indices can be written.
 */
class ResWriter {
    private:
    std::ofstream file_path;

    public:
    explicit ResWriter(const std::string& filename = "");
    ~ResWriter();

    ResWriter(ResWriter&& other) noexcept;
    ResWriter& operator=(ResWriter&& other) noexcept;

    ResWriter(const ResWriter&) = delete;
    ResWriter& operator=(const ResWriter&) = delete;

    void open(const std::string& filename);
    void close();

    void add_loadcase(int id, WriterStepType step_type = WriterStepType::Static);

    void write_field(const model::Field& field,
                     const std::string& field_name,
                     const model::ModelData* model_data = nullptr,
                     Precision frame_value = std::numeric_limits<Precision>::quiet_NaN());
};

} // namespace writer
} // namespace io
} // namespace fem
