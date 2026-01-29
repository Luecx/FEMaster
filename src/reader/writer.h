#pragma once

#include "../core/logging.h"
#include "../core/core.h"

#include <Eigen/Dense>
#include <fstream>
#include <iostream>

namespace fem {
namespace model {
struct Field;
}
namespace reader {

// Import DynamicMatrix from the outer fem namespace
using fem::DynamicMatrix;

/**
 * @brief Simple result file writer for FEMaster .res format.
 */
class Writer {
    private:
    std::ofstream file_path;

    public:
    // Constructor
    explicit Writer(const std::string& filename);

    // Destructor
    ~Writer();

    // Move constructor
    Writer(Writer&& other) noexcept
        : file_path(std::move(other.file_path)) {}

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

    /**
     * @brief Write any Eigen matrix as a FIELD block.
     *
     * If index_cols == 0 (default), the legacy format is used:
     *   FIELD, NAME=..., COLS=..., ROWS=...
     *
     * If index_cols > 0, the new indexed format is used:
     *   FIELD, NAME=..., INDEX_COLS=..., VALUE_COLS=..., ROWS=...
     *
     * where VALUE_COLS = matrix.cols() - index_cols.
     */
    void write_eigen_matrix(const DynamicMatrix& matrix,
                            const std::string& field_name,
                            int index_cols = 0);

    /**
     * @brief Write a model::Field as a FIELD block.
     *
     * Mirrors write_eigen_matrix but reads values directly from the field.
     */
    void write_field(const model::Field& field,
                     const std::string& field_name,
                     int index_cols = 0);
};

} // namespace reader
} // namespace fem
