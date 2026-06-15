// Created by Luecx on 06.09.2023.
//

#include "writer.h"

#include "../core/core.h"
#include "../core/logging.h"
#include "../data/field.h"

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
 * @brief Infers the result domain for matrix output when index columns are
 *        present but no explicit domain was requested.
 */
model::FieldDomain result_matrix_domain(model::FieldDomain requested,
                                         int index_cols) {
    if (requested != model::FieldDomain::UNKNOWN) {
        return requested;
    }

    return index_cols > 0 ? model::FieldDomain::ELEMENT_NODAL
                          : model::FieldDomain::UNKNOWN;
}

/**
 * @brief Writes the common field header.
 *
 * Format:
 *
 * @code
 * FIELD, NAME=<name>, TYPE=<domain>, COLS=<values>, ROWS=<rows>
 * FIELD, NAME=<name>, TYPE=<domain>, INDEX_COLS=<indices>, VALUE_COLS=<values>, ROWS=<rows>
 * @endcode
 */
void write_field_header(std::ostream& os,
                        const std::string& field_name,
                        model::FieldDomain domain,
                        Index rows,
                        Index index_cols,
                        Index value_cols) {
    os << "FIELD, NAME=" << field_name
       << ", TYPE=" << field_domain_to_string(domain);

    if (index_cols == 0) {
        os << ", COLS=" << value_cols;
    } else {
        os << ", INDEX_COLS=" << index_cols
           << ", VALUE_COLS=" << value_cols;
    }

    os << ", ROWS=" << rows << '\n';

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
 * @brief Writes a field without leading index columns.
 */
void write_unindexed_field(std::ostream& os,
                           const model::Field& field) {
    for (Index row = 0; row < field.rows; ++row) {
        write_field_values(os, field, row);
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
    logging::error(file_path.is_open(),
                   "Cannot write field '", field_name,
                   "': file is not open.");

    const int rows = static_cast<int>(matrix.rows());
    const int cols = static_cast<int>(matrix.cols());

    logging::error(index_cols >= 0 && index_cols <= cols,
                   "Writer: field '", field_name,
                   "' has invalid index column count");

    const int value_cols = cols - index_cols;
    logging::error(index_cols == 0 || value_cols > 0,
                   "Writer: field '", field_name,
                   "' must have at least one value column");

    const auto domain = result_matrix_domain(field_type, index_cols);

    write_field_header(file_path,
                       field_name,
                       domain,
                       static_cast<Index>(rows),
                       static_cast<Index>(index_cols),
                       static_cast<Index>(value_cols));

    for (int row = 0; row < rows; ++row) {
        write_matrix_values(file_path,
                            matrix,
                            row);
    }

    write_field_footer(file_path);
}

void Writer::write_field(const model::Field& field,
                         const std::string& field_name,
                         const model::ModelData* model_data) {
    (void) model_data;

    logging::error(file_path.is_open(),
                   "Cannot write field '", field_name,
                   "': file is not open.");

    write_field_header(file_path,
                       field_name,
                       field.domain,
                       field.rows,
                       0,
                       field.components);

    write_unindexed_field(file_path,
                          field);

    write_field_footer(file_path);
}

} // namespace reader
} // namespace fem
