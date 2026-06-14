#pragma once

#include "../core/logging.h"
#include "../core/core.h"
#include "../data/field.h"
#include "../core/types_cls.h"

#include <Eigen/Dense>
#include <fstream>
#include <iostream>

namespace fem {
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
     * Dense matrix fallback for non-topological result blocks.
     */
    void write_eigen_matrix(const DynamicMatrix& matrix,
                            const std::string& field_name,
                            int index_cols = 0,
                            model::FieldDomain field_type = model::FieldDomain::UNKNOWN);

    /**
     * @brief Write a model::Field as a FIELD block.
     *
     * Writes a field with indices derived from its domain. NODE and ELEMENT
     * rows get one index, ELEMENT_NODAL and ELEMENT_IP rows get element/local
     * indices. Rows with all value components equal to zero are omitted.
     */
    void write_field(const model::Field& field,
                     const std::string& field_name,
                     const model::ModelData* model_data = nullptr);
};
} // namespace reader
} // namespace fem
