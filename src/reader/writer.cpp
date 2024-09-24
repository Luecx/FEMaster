//
// Created by Luecx on 06.09.2023.
//
#include "writer.h"
#include "../core/core.h"

fem::reader::Writer::Writer(const std::string& filename) {
    this->open(filename);
}
fem::reader::Writer::~Writer() {
    if (file_path.is_open()) {
        file_path.close();
    }
}
fem::reader::Writer& fem::reader::Writer::operator=(fem::reader::Writer&& other) noexcept {
    if (this != &other) {
        close();
        file_path = std::move(other.file_path);
    }
    return *this;
}
void fem::reader::Writer::open(const std::string& filename) {
    close();
    file_path.open(filename);
    logging::error(file_path.is_open(), "Failed to open file: ", filename);
}
void fem::reader::Writer::close() {
    if (file_path.is_open()) {
        file_path.close();
    }
}
void fem::reader::Writer::add_loadcase(int id) {
    file_path << "LC " << id << std::endl;
}
void fem::reader::Writer::write_eigen_matrix(const DynamicMatrix& matrix, const std::string& field_name) {
    file_path << "FIELD, NAME=" << field_name << ", COLS=" << matrix.cols() << ", ROWS=" << matrix.rows() << std::endl;
    file_path << matrix << std::endl;
    file_path << "END FIELD" << std::endl;
}
