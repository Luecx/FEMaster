#pragma once

#include "../core/core.h"
#include "../core/logging.h"
#include "../data/field.h"

#include <fstream>
#include <string>
#include <vector>

namespace fem {
namespace model {
struct ModelData;
struct ElementInterface;
}

namespace reader {

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

    void add_loadcase(int id);

    void write_model_data(const model::ModelData& model_data);
    void write_field(const model::Field& field,
                     const std::string& field_name,
                     const model::ModelData* model_data = nullptr);

    private:
    void write_nodes      (const model::ModelData& model_data);
    void write_elements   (const model::ModelData& model_data);
    void write_nodal_field(const model::Field& field,
                           const std::string& field_name);

    ID  frd_node_number   (ID internal_node_id) const;
    ID  frd_element_number(ID internal_element_id) const;
};

} // namespace reader
} // namespace fem