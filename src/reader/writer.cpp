//
// Created by Luecx on 06.09.2023.
//
#include "writer.h"
#include "../core/core.h"

fem::reader::Writer::Writer(const std::string& filename) {
    this->open(filename);
}
fem::reader::Writer::~Writer() {
    if (out_file.is_open()) {
        out_file.close();
    }
}
fem::reader::Writer& fem::reader::Writer::operator=(fem::reader::Writer&& other) noexcept {
    if (this != &other) {
        close();
        out_file = std::move(other.out_file);
    }
    return *this;
}
void fem::reader::Writer::open(const std::string& filename) {
    close();
    out_file.open(filename);
    logging::error(out_file.is_open(), "Failed to open file: ", filename);
}
void fem::reader::Writer::close() {
    if (out_file.is_open()) {
        out_file.close();
    }
}
void fem::reader::Writer::add_loadcase(int id) {
    out_file << "LC " << id << std::endl;
}
void fem::reader::Writer::write_eigen_matrix(const DynamicMatrix& matrix, const std::string& field_name) {
    out_file << "FIELD, NAME=" << field_name << ", COLS=" << matrix.cols() << std::endl;
    out_file << matrix << std::endl;
    out_file << "END FIELD" << std::endl;
}
