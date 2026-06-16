// Created by Luecx on 06.09.2023.
//
#include "writer.h"
#include "../core/core.h"
#include "../core/logging.h"
#include "../data/field.h"
#include "../model/element/element.h"
#include "../model/model_data.h"

#include <cstddef>
#include <iomanip> // std::setw

namespace fem {
namespace reader {

namespace {

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

model::FieldDomain result_matrix_domain(model::FieldDomain requested, int index_cols) {
    if (requested != model::FieldDomain::UNKNOWN) {
        return requested;
    }
    return index_cols > 0 ? model::FieldDomain::ELEMENT_NODAL : model::FieldDomain::UNKNOWN;
}

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

void write_field_values(std::ofstream& file_path,
                        const model::Field& field,
                        Index row) {
    for (Index j = 0; j < field.components; ++j) {
        file_path << std::setw(16) << field(row, j);
    }
    file_path << '\n';
}

void write_element_location_field(std::ofstream& file_path,
                                  const model::Field& field,
                                  const std::string& field_name,
                                  const model::ModelData& model_data,
                                  const model::Field& offsets,
                                  const char* offset_name) {
    logging::error(offsets.rows == static_cast<Index>(model_data.max_elems + 1),
                   "Writer: ", offset_name, " must have max_elems + 1 rows");
    logging::error(field.rows == static_cast<Index>(offsets(static_cast<Index>(model_data.max_elems), 0)),
                   "Writer: field '", field_name, "' row count does not match ", offset_name);

    write_indexed_field_header(file_path,
                               field_name,
                               field.domain,
                               2,
                               static_cast<int>(field.components),
                               static_cast<int>(field.rows));

    file_path.setf(std::ios::scientific, std::ios::floatfield);

    for (Index elem_row = 0; elem_row < static_cast<Index>(model_data.max_elems); ++elem_row) {
        const auto& element = model_data.elements[static_cast<std::size_t>(elem_row)];
        if (!element) {
            continue;
        }

        const Index begin = static_cast<Index>(offsets(elem_row, 0));
        const Index end = static_cast<Index>(offsets(elem_row + 1, 0));
        logging::error(begin <= end && end <= field.rows,
                       "Writer: invalid ", offset_name, " span for element ", elem_row);

        const ID elem_id = element->elem_id;
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

Writer::Writer(const std::string& filename) {
    if (!filename.empty()) {
        this->open(filename);
    }
}

Writer::~Writer() {
    if (file_path.is_open()) {
        file_path.close();
    }
}

Writer& Writer::operator=(Writer&& other) noexcept {
    if (this != &other) {
        close();
        file_path = std::move(other.file_path);
    }
    return *this;
}

void Writer::open(const std::string& filename) {
    close();
    file_path.open(filename);
    logging::error(file_path.is_open(), "Failed to open file: ", filename);
}

void Writer::close() {
    if (file_path.is_open()) {
        file_path.close();
    }
}

void Writer::add_loadcase(int id) {
    file_path << "LC " << id << std::endl;
}

void Writer::write_eigen_matrix(const DynamicMatrix& matrix,
                                const std::string& field_name,
                                int index_cols,
                                model::FieldDomain field_type) {
    if (!file_path.is_open()) {
        logging::error(false, "Cannot write field '", field_name,
                       "': file is not open.");
    }

    const int rows = static_cast<int>(matrix.rows());
    const int cols = static_cast<int>(matrix.cols());

    logging::error(index_cols >= 0 && index_cols <= cols,
                   "Invalid index_cols for field ", field_name,
                   ": ", index_cols, " (cols=", cols, ")");

    const auto domain = result_matrix_domain(field_type, index_cols);

    if (index_cols == 0) {
        write_dense_field_header(file_path, field_name, domain, cols, rows);
    } else {
        const int value_cols = cols - index_cols;
        logging::error(value_cols > 0,
                       "Field ", field_name,
                       " must have at least one value column (index_cols=",
                       index_cols, ", cols=", cols, ")");

        write_indexed_field_header(file_path, field_name, domain, index_cols, value_cols, rows);
    }

    file_path.setf(std::ios::scientific, std::ios::floatfield);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file_path << std::setw(16) << matrix(i, j);
        }
        file_path << '\n';
    }

    file_path << "END FIELD\n";
}

void Writer::write_field(const model::Field& field,
                         const std::string& field_name,
                         int index_cols) {
    if (!file_path.is_open()) {
        logging::error(false, "Cannot write field '", field_name,
                       "': file is not open.");
    }

    const int rows = static_cast<int>(field.rows);
    const int cols = static_cast<int>(field.components);

    logging::error(index_cols >= 0 && index_cols <= cols,
                   "Invalid index_cols for field ", field_name,
                   ": ", index_cols, " (cols=", cols, ")");

    if (index_cols == 0) {
        write_dense_field_header(file_path, field_name, field.domain, cols, rows);
    } else {
        const int value_cols = cols - index_cols;
        logging::error(value_cols > 0,
                       "Field ", field_name,
                       " must have at least one value column (index_cols=",
                       index_cols, ", cols=", cols, ")");

        write_indexed_field_header(file_path, field_name, field.domain, index_cols, value_cols, rows);
    }

    file_path.setf(std::ios::scientific, std::ios::floatfield);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            file_path << std::setw(16) << field(i, j);
        }
        file_path << '\n';
    }

    file_path << "END FIELD\n";
}

void Writer::write_field(const model::Field& field,
                         const std::string& field_name,
                         const model::ModelData* model_data) {
    if (!file_path.is_open()) {
        logging::error(false, "Cannot write field '", field_name,
                       "': file is not open.");
    }

    if (field.domain == model::FieldDomain::ELEMENT_NODAL) {
        logging::error(model_data != nullptr,
                       "Writer: ELEMENT_NODAL field '", field_name, "' requires model data");
        logging::error(model_data->element_nodal_offsets != nullptr,
                       "Writer: element nodal offsets are not initialized for field '", field_name, "'");
        write_element_location_field(file_path,
                                     field,
                                     field_name,
                                     *model_data,
                                     *model_data->element_nodal_offsets,
                                     "element nodal offsets");
        return;
    }

    if (field.domain == model::FieldDomain::ELEMENT_IP) {
        logging::error(model_data != nullptr,
                       "Writer: ELEMENT_IP field '", field_name, "' requires model data");
        logging::error(model_data->element_ip_offsets != nullptr,
                       "Writer: element IP offsets are not initialized for field '", field_name, "'");
        write_element_location_field(file_path,
                                     field,
                                     field_name,
                                     *model_data,
                                     *model_data->element_ip_offsets,
                                     "element IP offsets");
        return;
    }

    write_field(field, field_name, 0);
}

} // namespace reader
} // namespace fem
