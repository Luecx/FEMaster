/**
 * @file writer_res.h
 * @brief Declares the native FEMaster result writer.
 *
 * The RES writer stores model fields together with their semantic entity
 * identifiers. NODE and ELEMENT fields are written as indexed fields instead of
 * relying on the dense solver row number as an implicit identifier.
 *
 * Entities of the implicit default Instance retain their original part-local
 * identifiers. Entities of explicit Instances use the qualified form
 * `instance_name.local_id`.
 *
 * Element-location fields additionally retain their local node, integration
 * point or material-point index.
 *
 * @see ResWriter
 * @see model::ModelData
 *
 * @author Finn Eggers
 * @date 20.08.2026
 */

#pragma once

#include "../../core/core.h"
#include "../../core/logging.h"
#include "../../data/field.h"
#include "writer_step_type.h"

#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace fem {
namespace model {
struct ModelData;
}

namespace io {
namespace writer {

/**
 * @brief FEMaster native `.res` result writer.
 *
 * NODE and ELEMENT fields are written with explicit semantic identifiers.
 * Entities of the default Instance retain their original numeric input
 * identifier, while entities of explicit Instances are qualified as
 * `instance_name.local_id`.
 *
 * Dense-to-semantic node and element identifiers are constructed once from the
 * compiled ModelData and cached by the writer for subsequent result fields.
 *
 * ELEMENT_NODAL, ELEMENT_IP and ELEMENT_MP fields use the semantic element
 * identifier together with their existing local-location index.
 */
class ResWriter {
    private:
    std::ofstream file_path;

    // Dense assembly ids -> semantic RES identifiers
    std::vector<std::string> node_ids;
    std::vector<std::string> element_ids;

    bool ids_initialized = false;

    public:
    explicit ResWriter(const std::string& filename = "");
    ~ResWriter();

    ResWriter(ResWriter&& other) noexcept;
    ResWriter& operator=(ResWriter&& other) noexcept;

    ResWriter(const ResWriter&) = delete;
    ResWriter& operator=(const ResWriter&) = delete;

    // File lifetime
    void open(const std::string& filename);
    void close();

    // Result organization
    void add_loadcase(int id, WriterStepType step_type = WriterStepType::Static);
    void write_pretension_sections(const model::ModelData& model_data);

    // Field output
    void write_field(const model::Field& field,
                     const std::string& field_name,
                     const model::ModelData* model_data = nullptr,
                     Precision frame_value = std::numeric_limits<Precision>::quiet_NaN());

    private:
    // Build semantic node and element identifiers from compiled Instance maps
    void initialize_ids(const model::ModelData& model_data);
};

} // namespace writer
} // namespace io
} // namespace fem
