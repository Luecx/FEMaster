#pragma once

#include "../core/logging.h"
#include "../core/core.h"

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
    Writer(const std::string& filename);

    // Destructor
    ~Writer();

    // Move constructor
    Writer(Writer&& other) noexcept
            : out_file(std::move(other.out_file)) {
    }

    // Move assignment operator
    Writer& operator=(Writer&& other) noexcept;

    // Delete copy constructor and copy assignment operator to prevent copying
    Writer(const Writer&) = delete;
    Writer& operator=(const Writer&) = delete;

    // Function to open the file
    void open(const std::string& filename);

    // Function to close the file
    void close();

    // Function to add a loadcase string
    void add_loadcase(int id);

    // Function to write any eigen matrix
    void write_eigen_matrix(const DynamicMatrix& matrix,
                            const std::string& field_name);
};

}    // namespace reader
}    // namespace fem
