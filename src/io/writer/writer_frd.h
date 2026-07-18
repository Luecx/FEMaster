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
struct ElementInterface;
}

namespace io {
namespace writer {

/**
 * @brief Minimal ASCII CalculiX/CGX .frd writer.
 *
 * Scope:
 * - writes node coordinates
 * - writes supported element connectivities
 * - writes NODE fields as nodal result blocks
 *
 * Non-nodal fields are ignored intentionally.
 */
class FrdWriter {
    private:
    std::ofstream file_path;

    bool model_data_written = false;

    int current_step         = 1;
    int current_result_block = 0;
    int current_global_frame = 0;
    int last_frame_step      = -1;
    int last_frame_id        = -1;

    WriterStepType current_step_type = WriterStepType::Static;

    bool remap_zero_node_ids    = false;
    bool remap_zero_element_ids = false;

    std::vector<ID> frd_node_ids;

    public:
    explicit FrdWriter(const std::string& filename = "");
    ~FrdWriter();

    FrdWriter(const FrdWriter&) = delete;
    FrdWriter& operator=(const FrdWriter&) = delete;

    void open(const std::string& filename);
    void close();

    void add_loadcase(int id, WriterStepType step_type = WriterStepType::Static);

    void write_model_data(const model::ModelData& model_data);
    void write_field(const model::Field& field,
                     const std::string& field_name,
                     const model::ModelData* model_data = nullptr,
                     Precision frame_value = std::numeric_limits<Precision>::quiet_NaN());

    private:
    void write_nodes      (const model::ModelData& model_data);
    void write_elements   (const model::ModelData& model_data);
    void write_nodal_field(const model::Field& field,
                           const std::string& field_name,
                           Precision frame_value);

    ID  frd_node_number   (ID internal_node_id) const;
    ID  frd_element_number(ID internal_element_id) const;
};

} // namespace writer
} // namespace io
} // namespace fem
