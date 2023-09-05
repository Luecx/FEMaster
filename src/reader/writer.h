#pragma once

#include "../core/logging.h"

#include <Eigen/Dense>
#include <fstream>
#include <iostream>

namespace fem {
namespace reader {

class writer {
    private:
    std::ofstream out_file;

    public:
    // Constructor
    writer() {}

    // Destructor
    ~writer() {
        if (out_file.is_open()) {
            out_file.close();
        }
    }

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
    template<typename T>
    void write_eigen_matrix(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix,
                            const std::string&                                      field_name) {
        out_file << "FIELD, NAME=" << field_name << ", COLS=" << matrix.cols() << std::endl;
        for (int i = 0; i < matrix.rows(); ++i) {
            for (int j = 0; j < matrix.cols(); ++j) {
                out_file << matrix(i, j) << " ";
            }
            out_file << std::endl;
        }
    }
};

}    // namespace reader
}    // namespace fem