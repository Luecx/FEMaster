#pragma once

#include "../core/logging.h"

#include <Eigen/Dense>
#include <fstream>
#include <iostream>

namespace fem {
namespace reader {

class Writer {
private:
    std::ofstream out_file;

public:
    // Constructor
    Writer(const std::string& filename) {
        this->open(filename);
    }

    // Destructor
    ~Writer() {
        if (out_file.is_open()) {
            out_file.close();
        }
    }

    // Move constructor
    Writer(Writer&& other) noexcept
            : out_file(std::move(other.out_file)) {
    }

    // Move assignment operator
    Writer& operator=(Writer&& other) noexcept {
        if (this != &other) {
            close();
            out_file = std::move(other.out_file);
        }
        return *this;
    }

    // Delete copy constructor and copy assignment operator to prevent copying
    Writer(const Writer&) = delete;
    Writer& operator=(const Writer&) = delete;

    // Function to open the file
    void open(const std::string& filename) {
        close();
        out_file.open(filename);
        logging::error(out_file.is_open(), "Failed to open file: ", filename);
    }

    // Function to close the file
    void close() {
        if (out_file.is_open()) {
            out_file.close();
        }
    }

    // Function to add a loadcase string
    void add_loadcase(int id) {
        out_file << "LC " << id << std::endl;
    }

    // Function to write any eigen matrix
    void write_eigen_matrix(const DynamicMatrix& matrix,
                            const std::string& field_name) {
        out_file << "FIELD, NAME=" << field_name << ", COLS=" << matrix.cols() << std::endl;
        out_file << matrix << std::endl;
        out_file << "END FIELD" << std::endl;
    }
};

}    // namespace reader
}    // namespace fem
