// Created by Luecx on 06.09.2023.
//
#include "writer.h"
#include "../core/core.h"
#include "../core/logging.h"
#include "../data/field.h"

#include <iomanip> // std::setw

namespace fem {
namespace reader {

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
                                int index_cols) {
    if (!file_path.is_open()) {
        logging::error(false, "Cannot write field '", field_name,
                       "': file is not open.");
    }

    const int rows = static_cast<int>(matrix.rows());
    const int cols = static_cast<int>(matrix.cols());

    logging::error(index_cols >= 0 && index_cols <= cols,
                   "Invalid index_cols for field ", field_name,
                   ": ", index_cols, " (cols=", cols, ")");

    if (index_cols == 0) {
        // Legacy behavior: no explicit index/value split
        file_path << "FIELD, NAME=" << field_name
                  << ", COLS=" << cols
                  << ", ROWS=" << rows << "\n";
    } else {
        const int value_cols = cols - index_cols;
        logging::error(value_cols > 0,
                       "Field ", field_name,
                       " must have at least one value column (index_cols=",
                       index_cols, ", cols=", cols, ")");

        file_path << "FIELD, NAME=" << field_name
                  << ", INDEX_COLS=" << index_cols
                  << ", VALUE_COLS=" << value_cols
                  << ", ROWS=" << rows << "\n";
    }

    // Dump matrix data
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
        file_path << "FIELD, NAME=" << field_name
                  << ", COLS=" << cols
                  << ", ROWS=" << rows << "\n";
    } else {
        const int value_cols = cols - index_cols;
        logging::error(value_cols > 0,
                       "Field ", field_name,
                       " must have at least one value column (index_cols=",
                       index_cols, ", cols=", cols, ")");

        file_path << "FIELD, NAME=" << field_name
                  << ", INDEX_COLS=" << index_cols
                  << ", VALUE_COLS=" << value_cols
                  << ", ROWS=" << rows << "\n";
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

} // namespace reader
} // namespace fem
