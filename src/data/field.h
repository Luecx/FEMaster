/**
 * @file field.h
 * @brief Defines owned dense numerical storage and model-domain field metadata.
 *
 * The model-data subsystem stores field values in contiguous row-major order.
 * FieldMatrix manages dimensions, scalar storage and bounds-checked access.
 * Field adds a name, entity domain, row/component counts and componentwise
 * arithmetic. Model compilation and result recovery determine row meaning,
 * physical units, tensor conventions and coordinate bases.
 *
 * FieldDomain identifies the association of each row; it does not perform
 * identifier remapping. Field is independent of fem::Namable and owns mutable
 * metadata together with its numerical storage.
 *
 * @see FieldMatrix
 * @see Field
 * @see FieldDomain
 * @see ModelData
 *
 * @author Finn Eggers
 * @date 05.10.2026
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

/**
 * @brief Identifies the entity layout represented by numerical field rows.
 *
 * Domains distinguish nodes, elements and element-local nodes, integration points
 * or material points. Dense row ordering and global offsets are assigned by model
 * compilation; this tag alone does not identify the referenced entity or basis.
 */
enum class FieldDomain : std::uint8_t {
    // Unspecified row association.
    UNKNOWN,
    // One row per model node.
    NODE,
    // One row per model element.
    ELEMENT,
    // Rows for element-local nodes in the compiled offset layout.
    ELEMENT_NODAL,
    // Rows for element integration points in the compiled offset layout.
    ELEMENT_IP,
    // Rows for material points in the compiled offset layout.
    ELEMENT_MP
};

/**
 * @brief Owns a dense row-major matrix of Precision scalar values.
 *
 * Storage uses offset(row, col) = row * cols + col. Sized construction allocates
 * zero-initialized values; the default object is an empty matrix. Dimensions use
 * the unsigned Index type, and callers must ensure their product fits storage
 * size and available memory. Indexing checks both upper bounds through logging.
 *
 * Copies own independent scalar vectors. data() exposes non-owning contiguous
 * pointers whose lifetime follows the storage; assignment or object destruction
 * can invalidate them. Fill operations modify only values, and finite-value
 * queries distinguish empty/all-nonfinite storage from at least one finite value.
 * Physical interpretation, entity layout and tensor operations belong to Field
 * and its consumers rather than to this storage abstraction.
 */
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

    // Non-owning access to contiguous row-major values, valid while storage is retained.
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

    // Persistent dimensions and owned contiguous row-major scalar storage.
    Index                  rows_{};
    Index                  cols_{};
    std::vector<Precision> data_{};
};

/**
 * @brief Owns named numerical values associated with a compiled model domain.
 *
 * Each row contains components scalar values in contiguous row-major order.
 * Sized construction requires positive row and component counts and initializes
 * values to zero. Default construction gives an empty field with UNKNOWN domain.
 * The public rows/components metadata must remain consistent with values; public
 * access does not enforce that invariant after construction. Copies deep-copy
 * numerical storage, while Ptr supplies optional shared ownership of a Field.
 *
 * Domain and dimensions must agree for componentwise field arithmetic; names,
 * units, component meaning and coordinate bases are not compared. Consumers own
 * these physical conventions and row remapping. Multiplication/division are
 * componentwise operations, not tensor products or matrix solves. Operations
 * modify this field in place and retain its metadata. Division checks exact zero
 * denominators; field division can update earlier entries before a later failure.
 *
 * Scalar access requires one component. row_vec3()/row_vec6() copy the leading
 * components without reordering or changing basis. Despite its name, is_nan()
 * reports all nonfinite values, including infinities. Raw data pointers expose
 * storage without ownership and require the Field to remain alive.
 */
struct Field {
    // Optional shared ownership of the complete field and its owned values.
    using Ptr = std::shared_ptr<Field>;

    // Mutable identifier and semantic association of rows with model entities.
    std::string name{};
    FieldDomain domain{FieldDomain::UNKNOWN};

    // Public dimensions; callers must preserve consistency with values dimensions.
    Index rows{};
    Index components{};

    // Owned values with offset(row, component) = row * components + component.
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

    // Detect NaN or infinity at a checked row/component index.
    [[nodiscard]] bool is_nan(Index row, Index component) const;

    // Validates that every stored component is finite and identifies invalid
    // rows and components through the supplied diagnostic label.
    void check_finite(const std::string& label) const;

    // Returns the first three components of a row
    [[nodiscard]] Vec3 row_vec3(Index row) const;

    // Returns the first six components of a row
    [[nodiscard]] Vec6 row_vec6(Index row) const;

    // Applies scalar compound operations
    Field& operator+=(Precision scalar);
    Field& operator-=(Precision scalar);
    Field& operator*=(Precision scalar);
    Field& operator/=(Precision scalar);

    // Apply componentwise operations after matching domains and dimensions.
    // Names, physical units and coordinate bases remain caller responsibilities.
    Field& operator+=(const Field& other);
    Field& operator-=(const Field& other);
    Field& operator*=(const Field& other);
    Field& operator/=(const Field& other);

private:
    // Validates domain and dimensions for an element-wise operation
    void validate_compatible(const Field& other, const char* operation) const;
};

// Legacy field spelling; the alias does not impose a NODE domain.
using NodeData = Field;

} // namespace model
} // namespace fem
