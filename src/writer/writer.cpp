// Created by Luecx on 06.09.2023.
//

#include "writer.h"

#include "../core/core.h"
#include "../core/logging.h"
#include "../data/field.h"
#include "../model/element/element.h"
#include "../model/model_data.h"

#include <iomanip>
#include <ostream>

namespace fem {
namespace reader {

namespace {

/**
 * @brief Column width used for all numeric field output.
 *
 * The writer intentionally uses a fixed-width text format. This makes the
 * result file easier to inspect manually and keeps columns aligned for
 * post-processing scripts that rely on whitespace-separated values.
 */
constexpr int FIELD_WIDTH = 16;

/**
 * @brief Converts a field domain enum to the textual representation used in
 *        the result file header.
 *
 * The output string is part of the external result format and should therefore
 * remain stable unless the result reader is updated accordingly.
 */
const char* field_domain_to_string(model::FieldDomain domain) {
    switch (domain) {
        case model::FieldDomain::UNKNOWN:
            return "UNKNOWN";

        case model::FieldDomain::NODE:
            return "NODE";

        case model::FieldDomain::ELEMENT:
            return "ELEMENT";

        case model::FieldDomain::ELEMENT_NODAL:
            return "ELEMENT_NODAL";

        case model::FieldDomain::ELEMENT_IP:
            return "ELEMENT_IP";
    }

    return "UNKNOWN";
}

/**
 * @brief Returns true if a matrix row contains at least one non-zero value.
 *
 * Empty rows are skipped in the result file. This keeps sparse result fields
 * compact, for example when only a subset of nodes or elements carries data.
 */
bool matrix_row_has_values(const DynamicMatrix& matrix,
                           int row) {
    for (int col = 0; col < matrix.cols(); ++col) {
        if (matrix(row, col) != Precision(0)) {
            return true;
        }
    }

    return false;
}

/**
 * @brief Returns true if a field row contains at least one non-zero component.
 *
 * This is used consistently for all field domains so that zero-valued rows are
 * omitted from the output file.
 */
bool field_row_has_values(const model::Field& field,
                          Index row) {
    for (Index component = 0; component < field.components; ++component) {
        if (field(row, component) != Precision(0)) {
            return true;
        }
    }

    return false;
}

/**
 * @brief Writes the common field header.
 *
 * Format:
 *
 * @code
 * FIELD, NAME=<name>, TYPE=<domain>, COLS=<number_of_components>
 * @endcode
 */
void write_field_header(std::ostream& os,
                        const std::string& field_name,
                        model::FieldDomain domain,
                        Index cols) {
    os << "FIELD, NAME=" << field_name
       << ", TYPE=" << field_domain_to_string(domain)
       << ", COLS=" << cols
       << '\n';

    os.setf(std::ios::scientific, std::ios::floatfield);
}

/**
 * @brief Writes the common field footer.
 */
void write_field_footer(std::ostream& os) {
    os << "END FIELD\n";
}

/**
 * @brief Writes all components of one dense matrix row.
 *
 * No row index is written here. The caller is responsible for writing any
 * leading index columns before calling this function.
 */
void write_matrix_values(std::ostream& os,
                         const DynamicMatrix& matrix,
                         int row) {
    for (int col = 0; col < matrix.cols(); ++col) {
        os << std::setw(FIELD_WIDTH) << matrix(row, col);
    }

    os << '\n';
}

/**
 * @brief Writes all components of one model field row.
 *
 * No node, element, integration-point, or local-node index is written here.
 * This function only serializes the actual field values.
 */
void write_field_values(std::ostream& os,
                        const model::Field& field,
                        Index row) {
    for (Index component = 0; component < field.components; ++component) {
        os << std::setw(FIELD_WIDTH) << field(row, component);
    }

    os << '\n';
}

/**
 * @brief Writes one indexed field row.
 *
 * Used for fields that have exactly one leading index column, for example
 * nodal fields or element fields.
 *
 * Format:
 *
 * @code
 * <id> <component_0> <component_1> ...
 * @endcode
 */
void write_indexed_field_row(std::ostream& os,
                             const model::Field& field,
                             Index row,
                             ID id) {
    if (!field_row_has_values(field, row)) {
        return;
    }

    os << std::setw(FIELD_WIDTH) << id;
    write_field_values(os, field, row);
}

/**
 * @brief Writes one element-local field row.
 *
 * Used for fields that belong to sublocations of an element, such as element
 * nodal values or integration point values.
 *
 * Format:
 *
 * @code
 * <element_id> <local_id> <component_0> <component_1> ...
 * @endcode
 *
 * For ELEMENT_NODAL fields, local_id is the local node index.
 * For ELEMENT_IP fields, local_id is the local integration point index.
 */
void write_subindexed_field_row(std::ostream& os,
                                const model::Field& field,
                                Index row,
                                ID element_id,
                                Index local_id) {
    if (!field_row_has_values(field, row)) {
        return;
    }

    os << std::setw(FIELD_WIDTH) << element_id;
    os << std::setw(FIELD_WIDTH) << local_id;

    write_field_values(os, field, row);
}

/**
 * @brief Writes a field without leading index columns.
 *
 * This is mainly used for fields with UNKNOWN domain, where no stable mapping
 * to nodes, elements, element nodes, or integration points is available.
 */
void write_unindexed_field(std::ostream& os,
                           const model::Field& field) {
    for (Index row = 0; row < field.rows; ++row) {
        if (!field_row_has_values(field, row)) {
            continue;
        }

        write_field_values(os, field, row);
    }
}

/**
 * @brief Writes a field with one leading integer index column.
 *
 * This currently assumes that row index and external entity id are identical.
 * If FEMaster allows sparse or non-contiguous node/element ids, this function
 * should be replaced by node-/element-aware variants using ModelData.
 */
void write_indexed_field(std::ostream& os,
                         const model::Field& field) {
    for (Index row = 0; row < field.rows; ++row) {
        write_indexed_field_row(os,
                                field,
                                row,
                                static_cast<ID>(row));
    }
}

/**
 * @brief Returns the correct element offset field for ELEMENT_NODAL or
 *        ELEMENT_IP data.
 *
 * ELEMENT_NODAL and ELEMENT_IP fields are stored as flattened arrays. The
 * offset field maps each element id to the start and end row of its local data.
 */
const model::Field& get_element_offsets(const model::Field& field,
                                        const model::ModelData& model_data,
                                        const std::string& field_name) {
    const auto& offsets = field.domain == model::FieldDomain::ELEMENT_NODAL
        ? model_data.element_nodal_offsets
        : model_data.element_ip_offsets;

    logging::error(offsets != nullptr,
                   "Writer: element offsets not initialized for field '",
                   field_name, "'");

    return *offsets;
}

/**
 * @brief Validates that a flattened element-local field matches its offset
 *        table.
 *
 * The final offset entry stores the total number of rows required by all
 * elements. Therefore, the field row count must match that value exactly.
 */
void validate_element_subfield(const model::Field& field,
                               const model::Field& element_offsets,
                               const model::ModelData& model_data,
                               const std::string& field_name) {
    const Index expected_rows = static_cast<Index>(
        element_offsets(static_cast<Index>(model_data.max_elems), 0)
    );

    logging::error(field.rows == expected_rows,
                   "Writer: field '", field_name,
                   "' row count does not match element offsets");
}

/**
 * @brief Writes an ELEMENT_NODAL or ELEMENT_IP field.
 *
 * The field is expected to be stored as a flattened array over all elements:
 *
 * @code
 * rows = [
 *     element_0 local rows,
 *     element_1 local rows,
 *     ...
 * ]
 * @endcode
 *
 * The corresponding offset field defines the row range for each element id.
 * Each written row receives two leading index columns:
 *
 * @code
 * <element_id> <local_id> <field components...>
 * @endcode
 */
void write_element_subindexed_field(std::ostream& os,
                                    const model::Field& field,
                                    const model::ModelData& model_data,
                                    const std::string& field_name) {
    const model::Field& element_offsets = get_element_offsets(field,
                                                              model_data,
                                                              field_name);

    validate_element_subfield(field,
                              element_offsets,
                              model_data,
                              field_name);

    for (Index elem_idx = 0; elem_idx < static_cast<Index>(model_data.max_elems); ++elem_idx) {
        const auto& element = model_data.elements[elem_idx];

        if (!element) {
            continue;
        }

        const ID    element_id = element->elem_id;
        const Index start      = static_cast<Index>(element_offsets(static_cast<Index>(element_id),     0));
        const Index end        = static_cast<Index>(element_offsets(static_cast<Index>(element_id) + 1, 0));

        for (Index row = start; row < end; ++row) {
            write_subindexed_field_row(os,
                                       field,
                                       row,
                                       element_id,
                                       row - start);
        }
    }
}

} // namespace


Writer::Writer(const std::string& filename) {
    if (!filename.empty()) {
        open(filename);
    }
}

Writer::~Writer() {
    close();
}

Writer& Writer::operator=(Writer&& other) noexcept {
    if (this == &other) {
        return *this;
    }

    close();

    file_path = std::move(other.file_path);

    return *this;
}

void Writer::open(const std::string& filename) {
    close();

    file_path.open(filename);

    logging::error(file_path.is_open(),
                   "Failed to open file: ", filename);
}

void Writer::close() {
    if (file_path.is_open()) {
        file_path.close();
    }
}

void Writer::add_loadcase(int id) {
    file_path << "LC " << id << '\n';
}

void Writer::write_eigen_matrix(const DynamicMatrix& matrix,
                                const std::string& field_name,
                                int index_cols,
                                model::FieldDomain field_type) {
    (void) index_cols;

    logging::error(file_path.is_open(),
                   "Cannot write field '", field_name,
                   "': file is not open.");

    const int rows = static_cast<int>(matrix.rows());
    const int cols = static_cast<int>(matrix.cols());

    write_field_header(file_path,
                       field_name,
                       field_type,
                       static_cast<Index>(cols));

    for (int row = 0; row < rows; ++row) {
        if (!matrix_row_has_values(matrix, row)) {
            continue;
        }

        write_matrix_values(file_path,
                            matrix,
                            row);
    }

    write_field_footer(file_path);
}

void Writer::write_field(const model::Field& field,
                         const std::string& field_name,
                         const model::ModelData* model_data) {
    logging::error(file_path.is_open(),
                   "Cannot write field '", field_name,
                   "': file is not open.");

    write_field_header(file_path,
                       field_name,
                       field.domain,
                       field.components);

    switch (field.domain) {
        case model::FieldDomain::NODE:
        case model::FieldDomain::ELEMENT:
            write_indexed_field(file_path,
                                field);
            break;

        case model::FieldDomain::ELEMENT_NODAL:
        case model::FieldDomain::ELEMENT_IP:
            logging::error(model_data != nullptr,
                           "Writer: field '", field_name,
                           "' requires ModelData for ELEMENT_NODAL/ELEMENT_IP output");

            write_element_subindexed_field(file_path,
                                           field,
                                           *model_data,
                                           field_name);
            break;

        case model::FieldDomain::UNKNOWN:
            write_unindexed_field(file_path,
                                  field);
            break;
    }

    write_field_footer(file_path);
}

} // namespace reader
} // namespace fem