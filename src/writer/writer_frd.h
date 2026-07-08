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
 * Non-nodal fields are ignored intentionally. ELEMENT_NODAL and ELEMENT_IP
 * should stay in the native .res output unless they are explicitly projected
 * or averaged to NODE fields before writing.
 */
class FrdWriter {
private:
    std::ofstream file_path;

    bool model_data_written = false;
    bool footer_written     = false;

    int current_loadcase = 0;
    int current_result_block = 0;

    bool remap_zero_node_ids    = false;
    bool remap_zero_element_ids = false;

    std::vector<ID> frd_node_ids;

public:
    explicit FrdWriter(const std::string& filename = "");
    ~FrdWriter();

    FrdWriter(FrdWriter&& other) noexcept;
    FrdWriter& operator=(FrdWriter&& other) noexcept;

    FrdWriter(const FrdWriter&) = delete;
    FrdWriter& operator=(const FrdWriter&) = delete;

    void open(const std::string& filename);
    void close();

    bool has_model_data() const;

    void write_model_data(const model::ModelData& model_data);

    void add_loadcase(int id);

    void write_field(const model::Field& field,
                     const std::string& field_name,
                     const model::ModelData* model_data = nullptr);

private:
    void write_header(const std::string& filename);
    void write_footer();

    void write_nodes(const model::ModelData& model_data);
    void write_elements(const model::ModelData& model_data);

    void write_nodal_field(const model::Field& field,
                           const std::string& field_name);

    bool supports_field(const model::Field& field,
                        const std::string& field_name) const;

    int frd_element_type(const model::ElementInterface& element) const;

    std::vector<ID> collect_frd_node_ids(const model::ModelData& model_data) const;

    std::size_t count_supported_elements(const model::ModelData& model_data) const;

    ID frd_node_number(ID internal_node_id) const;
    ID frd_element_number(ID internal_element_id) const;
};

} // namespace reader
} // namespace fem