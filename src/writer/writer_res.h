#pragma once

#include "../core/core.h"
#include "../core/logging.h"
#include "../data/field.h"

#include <fstream>
#include <string>

namespace fem {
namespace model {
struct ModelData;
}

namespace reader {

/**
 * @brief FEMaster .res writer.
 *
 * Writes all model::Field instances to the native FEMaster result format.
 * The field domain is taken from field.domain. ELEMENT_NODAL and ELEMENT_IP
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

    void add_loadcase(int id);

    void write_field(const model::Field& field,
                     const std::string& field_name,
                     const model::ModelData* model_data = nullptr);
};

} // namespace reader
} // namespace fem