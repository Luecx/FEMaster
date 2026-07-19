/**
 * @file field.h
 * @brief Declares lightweight dense storage for model fields.
 *
 * `Field` stores named scalar or vector-valued model data on nodes, elements,
 * integration points and other model domains. Values are stored contiguously
 * in row-major order as `rows x components`.
 *
 * Field semantics remain intentionally external. The container only manages
 * metadata, dense storage and basic element-wise operations.
 *
 * @see field.cpp
 * @see src/model/model_data.h
 */

#pragma once

#include "../core/types_eig.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace fem {
namespace model {

// Identifies the model entity associated with each field row
enum class FieldDomain : std::uint8_t {
    UNKNOWN,
    NODE,
    ELEMENT,
    ELEMENT_NODAL,
    ELEMENT_IP,
    ELEMENT_MP
};

// Dense row-major storage with index layout row * cols + col
class FieldMatrix {
public:
    // Constructs an empty matrix
    FieldMatrix() = default;

    // Constructs a zero-initialized matrix with the given dimensions
    FieldMatrix(Index rows, Index cols);

    // Returns the matrix dimensions
    [[nodiscard]] Index rows() const;
    [[nodiscard]] Index cols() const;

    // Returns the total number of stored scalar values
    [[nodiscard]] std::size_t size() const;

    // Returns a matrix entry by mutable or constant access
    Precision& operator()(Index row, Index col);
    Precision  operator()(Index row, Index col) const;

    // Returns the contiguous row-major storage
    [[nodiscard]] Precision*       data();
    [[nodiscard]] const Precision* data() const;

    // Initializes all stored values
    void set_zero();
    void set_ones();
    void fill_nan();

    // Returns whether at least one stored value is finite
    [[nodiscard]] bool has_any_finite() const;

private:
    // Converts a checked matrix index into a flat storage index
    [[nodiscard]] std::size_t offset(Index row, Index col) const;

    Index                  rows_{};
    Index                  cols_{};
    std::vector<Precision> data_{};
};

// Named dense field associated with a model domain
struct Field {
    using Ptr = std::shared_ptr<Field>;

    std::string name{};
    FieldDomain domain{FieldDomain::UNKNOWN};

    Index rows{};
    Index components{};

    FieldMatrix values{};

    // Constructs an empty field
    Field() = default;

    // Constructs an allocated field with the given metadata
    Field(std::string field_name,
          FieldDomain field_domain,
          Index       row_count,
          Index       component_count);

    // Returns a field component by mutable or constant access
    Precision& operator()(Index row, Index component);
    Precision  operator()(Index row, Index component) const;

    // Returns a scalar field value
    Precision& operator()(Index row);
    Precision  operator()(Index row) const;

    // Returns the contiguous row-major storage
    [[nodiscard]] Precision*       data();
    [[nodiscard]] const Precision* data() const;

    // Initializes all field values
    void set_zero();
    void set_ones();
    void fill_nan();

    // Returns whether at least one field value is finite
    [[nodiscard]] bool has_any_finite() const;

    // Returns whether a field component is not finite
    [[nodiscard]] bool is_nan(Index row, Index component) const;

    // Returns the first three components of a row
    [[nodiscard]] Vec3 row_vec3(Index row) const;

    // Returns the first six components of a row
    [[nodiscard]] Vec6 row_vec6(Index row) const;

    // Applies scalar compound operations
    Field& operator+=(Precision scalar);
    Field& operator-=(Precision scalar);
    Field& operator*=(Precision scalar);
    Field& operator/=(Precision scalar);

    // Applies element-wise compound operations
    Field& operator+=(const Field& other);
    Field& operator-=(const Field& other);
    Field& operator*=(const Field& other);
    Field& operator/=(const Field& other);

private:
    // Validates domain and dimensions for an element-wise operation
    void validate_compatible(const Field& other, const char* operation) const;
};

using NodeData = Field;

} // namespace model
} // namespace fem
