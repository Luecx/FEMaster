/**
 * @file writer_res.cpp
 * @brief Implements the native FEMaster result writer.
 *
 * NODE and ELEMENT fields are written as explicitly indexed fields. The writer
 * converts dense assembly identifiers to semantic result identifiers once from
 * the compiled Part/Instance topology and caches the resulting strings.
 *
 * The implicit default Instance preserves the original part-local identifier:
 *
 *     17
 *
 * Explicit Instances use a qualified identifier:
 *
 *     bolt.17
 *
 * The same element identifiers are reused for ELEMENT, ELEMENT_NODAL,
 * ELEMENT_IP and ELEMENT_MP fields.
 *
 * @see ResWriter
 * @see model::ModelData
 *
 * @author Finn Eggers
 * @date 20.08.2026
 */

#include "writer_res.h"

#include "../../model/element/element.h"
#include "../../model/model_data.h"

#include <cstddef>
#include <iomanip>
#include <utility>

namespace fem {
namespace io {
namespace writer {

namespace {

/**
 * Converts a field domain to its textual RES representation.
 */
const char* field_domain_to_string(model::FieldDomain domain) {
    switch (domain) {
        case model::FieldDomain::UNKNOWN:       return "UNKNOWN";
        case model::FieldDomain::NODE:          return "NODE";
        case model::FieldDomain::ELEMENT:       return "ELEMENT";
        case model::FieldDomain::ELEMENT_NODAL: return "ELEMENT_NODAL";
        case model::FieldDomain::ELEMENT_IP:    return "ELEMENT_IP";
        case model::FieldDomain::ELEMENT_MP:    return "ELEMENT_MP";
    }

    return "UNKNOWN";
}

/**
 * Writes the header of a dense RES field.
 */
void write_dense_field_header(std::ofstream& file_path,
                              const std::string& field_name,
                              model::FieldDomain domain,
                              int cols,
                              int rows) {
    file_path << "FIELD, NAME=" << field_name
              << ", TYPE=" << field_domain_to_string(domain)
              << ", COLS=" << cols
              << ", ROWS=" << rows << "\n";
}

/**
 * Writes the header of an explicitly indexed RES field.
 */
void write_indexed_field_header(std::ofstream& file_path,
                                const std::string& field_name,
                                model::FieldDomain domain,
                                int index_cols,
                                int value_cols,
                                int rows) {
    file_path << "FIELD, NAME=" << field_name
              << ", TYPE=" << field_domain_to_string(domain)
              << ", INDEX_COLS=" << index_cols
              << ", VALUE_COLS=" << value_cols
              << ", ROWS=" << rows << "\n";
}

/**
 * Writes all numerical components of one field row.
 */
void write_field_values(std::ofstream& file_path,
                        const model::Field& field,
                        Index row) {
    for (Index j = 0; j < field.components; ++j) {
        file_path << std::setw(16) << field(row, j);
    }

    file_path << '\n';
}

/**
 * Writes a field without explicit entity identifiers.
 *
 * This path is retained for UNKNOWN fields. All entity-based NODE and ELEMENT
 * domains are written through indexed output instead.
 */
void write_dense_field(std::ofstream& file_path,
                       const model::Field& field,
                       const std::string& field_name) {
    write_dense_field_header(
        file_path,
        field_name,
        field.domain,
        static_cast<int>(field.components),
        static_cast<int>(field.rows)
    );

    file_path.setf(std::ios::scientific, std::ios::floatfield);

    for (Index i = 0; i < field.rows; ++i) {
        write_field_values(file_path, field, i);
    }

    file_path << "END FIELD\n";
}

/**
 * Writes one NODE or ELEMENT field with an explicit semantic entity identifier.
 *
 * Field rows remain in dense assembly order. The corresponding semantic
 * identifier is supplied through `ids`, so no result interpretation depends on
 * the dense row number.
 */
void write_indexed_field(std::ofstream& file_path,
                         const model::Field& field,
                         const std::string& field_name,
                         const std::vector<std::string>& ids) {
    logging::error(ids.size() == static_cast<std::size_t>(field.rows),
        "ResWriter: id mapping does not match rows of field '", field_name, "'");

    write_indexed_field_header(
        file_path,
        field_name,
        field.domain,
        1,
        static_cast<int>(field.components),
        static_cast<int>(field.rows)
    );

    file_path.setf(std::ios::scientific, std::ios::floatfield);

    for (Index row = 0; row < field.rows; ++row) {
        file_path << std::setw(16)
                  << ids[static_cast<std::size_t>(row)];

        write_field_values(file_path, field, row);
    }

    file_path << "END FIELD\n";
}

/**
 * Writes an element-location field with semantic element identifiers.
 *
 * The first index identifies the element using its RES identifier. The second
 * index retains the existing local node, integration-point or material-point
 * enumeration within that element.
 */
void write_element_location_field(std::ofstream& file_path,
                                  const model::Field& field,
                                  const std::string& field_name,
                                  const model::ModelData& model_data,
                                  const std::vector<std::string>& element_ids,
                                  const model::Field& offsets,
                                  const char* offset_name) {
    const Index element_count = static_cast<Index>(model_data.elements.size());

    logging::error(offsets.rows == element_count + 1,
        "ResWriter: ", offset_name, " must have element_count + 1 rows");
    logging::error(field.rows == static_cast<Index>(offsets(element_count, 0)),
        "ResWriter: field '", field_name, "' row count does not match ", offset_name);
    logging::error(element_ids.size() == model_data.elements.size(),
        "ResWriter: element id mapping does not match compiled elements");

    write_indexed_field_header(
        file_path,
        field_name,
        field.domain,
        2,
        static_cast<int>(field.components),
        static_cast<int>(field.rows)
    );

    file_path.setf(std::ios::scientific, std::ios::floatfield);

    for (Index elem_row = 0; elem_row < element_count; ++elem_row) {
        const auto& element =
            model_data.elements[static_cast<std::size_t>(elem_row)];

        if (!element) {
            continue;
        }

        const Index begin = static_cast<Index>(offsets(elem_row, 0));
        const Index end   = static_cast<Index>(offsets(elem_row + 1, 0));

        logging::error(begin <= end && end <= field.rows,
            "ResWriter: invalid ", offset_name, " span for element ", elem_row);

        const std::string& elem_id =
            element_ids[static_cast<std::size_t>(elem_row)];

        for (Index local = 0; local < end - begin; ++local) {
            const Index row = begin + local;

            file_path << std::setw(16) << elem_id
                      << std::setw(16) << local;

            write_field_values(file_path, field, row);
        }
    }

    file_path << "END FIELD\n";
}

} // namespace

/**
 * Constructs a RES writer and optionally opens the target file.
 */
ResWriter::ResWriter(const std::string& filename) {
    if (!filename.empty()) {
        open(filename);
    }
}

/**
 * Closes the active result file.
 */
ResWriter::~ResWriter() {
    close();
}

/**
 * Moves an active writer together with its cached semantic entity identifiers.
 */
ResWriter::ResWriter(ResWriter&& other) noexcept
    : file_path      (std::move(other.file_path)),
      node_ids       (std::move(other.node_ids)),
      element_ids    (std::move(other.element_ids)),
      ids_initialized(other.ids_initialized) {
    other.ids_initialized = false;
}

/**
 * Moves an active writer together with its cached semantic entity identifiers.
 */
ResWriter& ResWriter::operator=(ResWriter&& other) noexcept {
    if (this != &other) {
        close();

        file_path       = std::move(other.file_path);
        node_ids        = std::move(other.node_ids);
        element_ids     = std::move(other.element_ids);
        ids_initialized = other.ids_initialized;

        other.ids_initialized = false;
    }

    return *this;
}

/**
 * Opens a new RES file and clears all writer-local identifier mappings.
 */
void ResWriter::open(const std::string& filename) {
    close();

    file_path.open(filename, std::ios::out | std::ios::trunc);

    logging::error(file_path.is_open(),
        "ResWriter: failed to open file: ", filename);

    node_ids.clear();
    element_ids.clear();
    ids_initialized = false;
}

/**
 * Closes the active RES file.
 */
void ResWriter::close() {
    if (file_path.is_open()) {
        file_path.close();
    }
}

/**
 * Starts a new result loadcase.
 */
void ResWriter::add_loadcase(int id, WriterStepType step_type) {
    (void) step_type;

    logging::error(file_path.is_open(),
        "ResWriter: cannot add loadcase: file is not open");

    file_path << "LC " << id << '\n';
}

/**
 * Builds the semantic RES identifiers for compiled nodes and elements.
 *
 * Nodes are resolved directly through `ModelData::node_mapping`, which maps a
 * dense assembly node to its owning Instance and original part-local node id.
 *
 * Elements are resolved by reversing each compiled Instance's existing
 * `element_ids` map from local to dense identifiers.
 *
 * The default Instance retains bare local identifiers for compatibility with
 * existing FEMaster result files. Explicit Instances are qualified as
 * `instance_name.local_id`.
 */
void ResWriter::initialize_ids(const model::ModelData& model_data) {
    logging::error(!ids_initialized,
        "ResWriter: result ids have already been initialized");

    // Build dense node -> semantic node identifier
    node_ids.resize(model_data.node_mapping.size());

    for (std::size_t row = 0; row < model_data.node_mapping.size(); ++row) {
        const auto& [instance, local_id] = model_data.node_mapping[row];

        logging::error(instance != nullptr,
            "ResWriter: dense node ", row, " has no source instance");

        node_ids[row] = instance->instance_id == 0
            ? std::to_string(local_id)
            : instance->name + "." + std::to_string(local_id);
    }

    // Reverse every Instance's local -> dense element map
    element_ids.assign(model_data.elements.size(), "");

    for (const auto& [instance_name, instance] : model_data.instances) {
        logging::error(instance != nullptr,
            "ResWriter: instance ", instance_name, " is null");

        for (const auto& [local_id, global_id] : instance->element_ids) {
            const std::size_t row = static_cast<std::size_t>(global_id);

            logging::error(row < element_ids.size(),
                "ResWriter: element ", local_id, " in instance ", instance_name,
                " maps outside compiled element storage");
            logging::error(element_ids[row].empty(),
                "ResWriter: dense element ", global_id,
                " has multiple semantic identifiers");

            element_ids[row] = instance->instance_id == 0
                ? std::to_string(local_id)
                : instance_name + "." + std::to_string(local_id);
        }
    }

    // Every compiled element must have a semantic source identifier
    for (std::size_t row = 0; row < model_data.elements.size(); ++row) {
        if (!model_data.elements[row]) {
            continue;
        }

        logging::error(!element_ids[row].empty(),
            "ResWriter: dense element ", row,
            " has no semantic element identifier");
    }

    ids_initialized = true;
}

/**
 * Writes one FEMaster result field.
 *
 * All entity-based field domains are explicitly indexed. The semantic node and
 * element mappings are initialized lazily on the first such write and retained
 * for subsequent fields.
 *
 * ELEMENT_NODAL, ELEMENT_IP and ELEMENT_MP fields still require ModelData on
 * every write because their dense row spans are defined by the corresponding
 * element-offset fields.
 */
void ResWriter::write_field(const model::Field& field,
                            const std::string& field_name,
                            const model::ModelData* model_data,
                            Precision frame_value) {
    (void) frame_value;

    logging::error(file_path.is_open(),
        "ResWriter: cannot write field '", field_name,
        "': file is not open");

    const bool indexed =
           field.domain == model::FieldDomain::NODE
        || field.domain == model::FieldDomain::ELEMENT
        || field.domain == model::FieldDomain::ELEMENT_NODAL
        || field.domain == model::FieldDomain::ELEMENT_IP
        || field.domain == model::FieldDomain::ELEMENT_MP;

    if (indexed && !ids_initialized) {
        logging::error(model_data != nullptr,
            "ResWriter: field '", field_name,
            "' requires model data to initialize result ids");

        initialize_ids(*model_data);
    }

    if (field.domain == model::FieldDomain::NODE) {
        write_indexed_field(file_path, field, field_name, node_ids);
        return;
    }

    if (field.domain == model::FieldDomain::ELEMENT) {
        write_indexed_field(file_path, field, field_name, element_ids);
        return;
    }

    if (field.domain == model::FieldDomain::ELEMENT_NODAL) {
        logging::error(model_data != nullptr,
            "ResWriter: ELEMENT_NODAL field '", field_name,
            "' requires model data");
        logging::error(model_data->element_nodal_offsets != nullptr,
            "ResWriter: element nodal offsets are not initialized for field '",
            field_name, "'");

        write_element_location_field(
            file_path,
            field,
            field_name,
            *model_data,
            element_ids,
            *model_data->element_nodal_offsets,
            "element nodal offsets"
        );
        return;
    }

    if (field.domain == model::FieldDomain::ELEMENT_IP) {
        logging::error(model_data != nullptr,
            "ResWriter: ELEMENT_IP field '", field_name,
            "' requires model data");
        logging::error(model_data->element_ip_offsets != nullptr,
            "ResWriter: element IP offsets are not initialized for field '",
            field_name, "'");

        write_element_location_field(
            file_path,
            field,
            field_name,
            *model_data,
            element_ids,
            *model_data->element_ip_offsets,
            "element IP offsets"
        );
        return;
    }

    if (field.domain == model::FieldDomain::ELEMENT_MP) {
        logging::error(model_data != nullptr,
            "ResWriter: ELEMENT_MP field '", field_name,
            "' requires model data");
        logging::error(model_data->element_mp_offsets != nullptr,
            "ResWriter: element MP offsets are not initialized for field '",
            field_name, "'");

        write_element_location_field(
            file_path,
            field,
            field_name,
            *model_data,
            element_ids,
            *model_data->element_mp_offsets,
            "element MP offsets"
        );
        return;
    }

    write_dense_field(file_path, field, field_name);
}

} // namespace writer
} // namespace io
} // namespace fem