#pragma once

#include "../core/logging.h"
#include "../core/core.h"
#include "../data/field.h"

#include <Eigen/Dense>
#include <fstream>
#include <iostream>

namespace fem {
namespace model {
struct ModelData;
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
     *   FIELD, NAME=..., TYPE=..., COLS=..., ROWS=...
     *
     * If index_cols > 0, the new indexed format is used:
     *   FIELD, NAME=..., TYPE=..., INDEX_COLS=..., VALUE_COLS=..., ROWS=...
     *
     * where VALUE_COLS = matrix.cols() - index_cols.
     * If no field_type is provided, indexed matrices are written as
     * ELEMENT_NODAL and dense matrices as UNKNOWN.
     */
    void write_eigen_matrix(const DynamicMatrix& matrix,
                            const std::string& field_name,
                            int index_cols = 0,
                            model::FieldDomain field_type = model::FieldDomain::UNKNOWN);

    /**
     * @brief Write a model::Field as a FIELD block.
     *
     * Mirrors write_eigen_matrix but reads values directly from the field.
     * The field type is taken from field.domain.
     */
    void write_field(const model::Field& field,
                     const std::string& field_name,
                     int index_cols = 0);

    /// Uses model topology to add element/local index columns for ELEMENT_NODAL and ELEMENT_IP fields.
    void write_field(const model::Field& field,
                     const std::string& field_name,
                     const model::ModelData* model_data);
};
} // namespace reader
} // namespace fem
