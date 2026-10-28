/**
 * @file field.cpp
 * @brief Implements dense field storage, validation and componentwise arithmetic.
 *
 * The model-data implementation allocates row-major scalar vectors, converts
 * checked row/component indices into storage offsets and exposes non-owning data
 * access. Field arithmetic checks domain and dimensions when combining fields.
 * Compound operators modify their left operand, while binary operators return
 * independent results retaining the metadata of the left operand.
 *
 * Finite-value diagnostics and zero-divisor checks report numerical failures
 * through logging. Units, tensor component ordering, coordinate bases and dense
 * entity-row mapping remain responsibilities of model and result consumers.
 *
 * @see FieldMatrix
 * @see Field
 * @see FieldDomain
 *
 * @author Finn Eggers
 * @date 05.10.2026
 */

#include "field.h"

#include "../core/logging.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace fem {
namespace model {

// FieldMatrix storage and indexing

/**
 * @brief Allocates zero-initialized contiguous matrix storage.
 *
 * The matrix owns `rows * cols` Precision values in row-major order. Zero-sized
 * dimensions are allowed; size multiplication is not checked for overflow.
 *
 * @param rows Number of rows in the matrix.
 * @param cols Number of components stored per row.
 */
FieldMatrix::FieldMatrix(Index rows, Index cols)
    : rows_(rows),
      cols_(cols) {
    // Convert dimensions to the vector storage type before multiplication.
    const auto row_count = static_cast<std::size_t>(rows);
    const auto col_count = static_cast<std::size_t>(cols);

    data_.resize(row_count * col_count);
}

/**
 * @brief Returns the number of matrix rows.
 *
 * The returned dimension is stored independently of the scalar vector.
 *
 * @return The configured row count.
 */
Index FieldMatrix::rows() const {
    return rows_;
}

/**
 * @brief Returns the number of matrix columns.
 *
 * Columns represent consecutive scalar entries within each row.
 *
 * @return The configured column count.
 */
Index FieldMatrix::cols() const {
    return cols_;
}

/**
 * @brief Returns the number of allocated scalar values.
 *
 * This is the actual storage size, not a recomputed dimension product.
 *
 * @return The number of stored Precision values.
 */
std::size_t FieldMatrix::size() const {
    return data_.size();
}

/**
 * @brief Returns mutable access to a checked matrix entry.
 *
 * The row and column are mapped to contiguous row-major storage after bounds
 * validation in `offset()`.
 *
 * @param row Zero-based row index.
 * @param col Zero-based column index.
 * @return Reference to the selected scalar.
 */
Precision& FieldMatrix::operator()(Index row, Index col) {
    return data_[offset(row, col)];
}

/**
 * @brief Reads a checked matrix entry.
 *
 * Uses the same bounds-checked row-major mapping as mutable access.
 *
 * @param row Zero-based row index.
 * @param col Zero-based column index.
 * @return Copy of the selected scalar.
 */
Precision FieldMatrix::operator()(Index row, Index col) const {
    return data_[offset(row, col)];
}

/**
 * @brief Exposes the mutable contiguous row-major storage.
 *
 * The pointer is non-owning and may be invalidated by storage replacement.
 *
 * @return Pointer to the first scalar, or the empty vector data pointer.
 */
Precision* FieldMatrix::data() {
    return data_.data();
}

/**
 * @brief Exposes the read-only contiguous row-major storage.
 *
 * No copy is made; the pointer is valid only while the storage remains alive.
 *
 * @return Const pointer to the scalar data.
 */
const Precision* FieldMatrix::data() const {
    return data_.data();
}

/**
 * @brief Sets every stored scalar to zero.
 *
 * Dimensions and storage capacity remain unchanged.
 */
void FieldMatrix::set_zero() {
    std::fill(data_.begin(), data_.end(), Precision(0));
}

/**
 * @brief Sets every stored scalar to one.
 *
 * Dimensions and storage capacity remain unchanged.
 */
void FieldMatrix::set_ones() {
    std::fill(data_.begin(), data_.end(), Precision(1));
}

/**
 * @brief Marks all matrix entries as unavailable.
 *
 * Uses the quiet NaN representation of Precision without changing the shape.
 */
void FieldMatrix::fill_nan() {
    const Precision nan = std::numeric_limits<Precision>::quiet_NaN();

    std::fill(data_.begin(), data_.end(), nan);
}

/**
 * @brief Tests whether at least one stored scalar is finite.
 *
 * The scan stops at the first finite value. An empty matrix or one consisting
 * only of NaN and infinity returns false.
 *
 * @return True if any scalar satisfies `std::isfinite`.
 */
bool FieldMatrix::has_any_finite() const {
    // An early success avoids scanning the remaining storage.
    for (const Precision value : data_) {
        if (std::isfinite(static_cast<double>(value))) {
            return true;
        }
    }

    return false;
}

/**
 * @brief Checks indices and computes the row-major storage offset.
 *
 * The flat offset is `row * cols + col`. Both indices must be within the
 * configured matrix dimensions before storage is accessed.
 *
 * @param row Zero-based row index.
 * @param col Zero-based column index.
 * @return Position of the selected scalar in the contiguous vector.
 */
std::size_t FieldMatrix::offset(Index row, Index col) const {
    // Reject invalid coordinates before computing the flattened index.
    logging::error(row < rows_ && col < cols_,
        "FieldMatrix: index (", row, ", ", col, ") is outside ", rows_, "x", cols_);

    return static_cast<std::size_t>(row) *
           static_cast<std::size_t>(cols_) +
           static_cast<std::size_t>(col);
}

// Field metadata, access and arithmetic

/**
 * @brief Constructs a named, zero-initialized field in a model domain.
 *
 * The field owns a dense matrix of the requested shape. Both dimensions must
 * be positive. Values are not converted between units or coordinate bases.
 *
 * @param field_name Name used to identify the field.
 * @param field_domain Domain associated with its rows.
 * @param row_count Number of field rows.
 * @param component_count Number of scalar components per row.
 */
Field::Field(std::string field_name,
             FieldDomain field_domain,
             Index       row_count,
             Index       component_count)
    : name      (std::move(field_name)),
      domain    (field_domain),
      rows      (row_count),
      components(component_count),
      values    (rows, components) {
    // Sized fields require a nonempty rectangular layout.
    logging::error(rows > 0,
        "Field '", name, "': rows must be positive");
    logging::error(components > 0,
        "Field '", name, "': components must be positive");
}

/**
 * @brief Creates an independent zero-filled field with the same layout.
 *
 * Copies the domain, row count and component count but uses the requested name.
 * The source field and its values are not modified.
 *
 * @param field_name Name assigned to the new field.
 * @return New field with zero-initialized storage.
 */
Field Field::zeros_like(std::string field_name) const {
    return Field(std::move(field_name), domain, rows, components);
}

/**
 * @brief Creates an independent NaN-filled field with the same layout.
 *
 * Preserves the source domain and dimensions. Nonfinite entries can be used
 * to mark values that are not provided for a result field.
 *
 * @param field_name Name assigned to the new field.
 * @return New field with all components set to quiet NaN.
 */
Field Field::nan_like(std::string field_name) const {
    // Reuse the matching layout before marking every component unavailable.
    Field result = zeros_like(std::move(field_name));
    result.fill_nan();
    return result;
}

/**
 * @brief Returns mutable access to a nodal or element field component.
 *
 * FieldMatrix performs the bounds check; this accessor does not reinterpret
 * the domain or component ordering.
 *
 * @param row Zero-based field row.
 * @param component Zero-based component.
 * @return Reference to the selected scalar.
 */
Precision& Field::operator()(Index row, Index component) {
    return values(row, component);
}

/**
 * @brief Reads one checked field component.
 *
 * The stored component is returned without unit or basis conversion.
 *
 * @param row Zero-based field row.
 * @param component Zero-based component.
 * @return Copy of the selected scalar.
 */
Precision Field::operator()(Index row, Index component) const {
    return values(row, component);
}

/**
 * @brief Returns mutable access to a scalar field row.
 *
 * This overload is only defined for fields with exactly one component.
 *
 * @param row Zero-based field row.
 * @return Reference to the row value.
 */
Precision& Field::operator()(Index row) {
    logging::error(components == 1,
        "Field '", name, "': scalar access requires exactly one component");

    return values(row, 0);
}

/**
 * @brief Reads the value of a scalar field row.
 *
 * Rejects multi-component fields so the requested component is unambiguous.
 *
 * @param row Zero-based field row.
 * @return Copy of the row value.
 */
Precision Field::operator()(Index row) const {
    logging::error(components == 1,
        "Field '", name, "': scalar access requires exactly one component");

    return values(row, 0);
}

/**
 * @brief Returns non-owning mutable access to the field storage.
 *
 * Entries are laid out row by row. The caller must preserve field lifetime.
 *
 * @return Pointer to the contiguous scalar data.
 */
Precision* Field::data() {
    return values.data();
}

/**
 * @brief Returns non-owning read-only access to the field storage.
 *
 * The returned pointer refers directly to the owned FieldMatrix values.
 *
 * @return Const pointer to the contiguous scalar data.
 */
const Precision* Field::data() const {
    return values.data();
}

/**
 * @brief Sets all field values to zero.
 *
 * The field name, domain and dimensions are preserved.
 */
void Field::set_zero() {
    values.set_zero();
}

/**
 * @brief Sets all field values to one.
 *
 * The field name, domain and dimensions are preserved.
 */
void Field::set_ones() {
    values.set_ones();
}

/**
 * @brief Sets all field values to quiet NaN.
 *
 * The field name, domain and dimensions are preserved.
 */
void Field::fill_nan() {
    values.fill_nan();
}

/**
 * @brief Checks whether the field contains at least one finite value.
 *
 * Delegates to the owned FieldMatrix without requiring all entries to be finite.
 *
 * @return True if any stored component is finite.
 */
bool Field::has_any_finite() const {
    return values.has_any_finite();
}

/**
 * @brief Tests a component for NaN or infinity.
 *
 * Despite its name, the method treats any nonfinite scalar as invalid.
 * Index validation is performed by the underlying matrix.
 *
 * @param row Zero-based field row.
 * @param component Zero-based component.
 * @return True if the selected value is not finite.
 */
bool Field::is_nan(Index row, Index component) const {
    const Precision value = values(row, component);

    return !std::isfinite(static_cast<double>(value));
}

/**
 * @brief Validates every field entry and reports its location on failure.
 *
 * The diagnostic includes the supplied label, row and component index.
 * Only values are checked; physical consistency and metadata/storage
 * invariants are outside the scope of this method.
 *
 * @param label Description of the quantity being validated.
 */
void Field::check_finite(const std::string& label) const {
    // Check each component so diagnostics identify its precise location.
    for (Index row = 0; row < rows; ++row) {
        for (Index component = 0; component < components; ++component) {
            logging::error(std::isfinite(static_cast<double>(values(row, component))),
                label, " row ", row, " has invalid value at component ", component);
        }
    }
}

/**
 * @brief Copies the leading three scalar components into a Vec3.
 *
 * At least three components are required. Additional components are ignored;
 * the stored component order and coordinate basis remain unchanged.
 *
 * @param row Zero-based field row.
 * @return Independent three-component vector.
 */
Vec3 Field::row_vec3(Index row) const {
    logging::error(components >= 3,
        "Field '", name, "': row_vec3 requires at least three components");

    return Vec3(
        values(row, 0),
        values(row, 1),
        values(row, 2)
    );
}

/**
 * @brief Copies the leading six scalar components into a Vec6.
 *
 * At least six components are required. Additional components are ignored;
 * the stored component order and coordinate basis remain unchanged.
 *
 * @param row Zero-based field row.
 * @return Independent six-component vector.
 */
Vec6 Field::row_vec6(Index row) const {
    logging::error(components >= 6,
        "Field '", name, "': row_vec6 requires at least six components");

    return Vec6(
        values(row, 0),
        values(row, 1),
        values(row, 2),
        values(row, 3),
        values(row, 4),
        values(row, 5)
    );
}

/**
 * @brief Adds a scalar to every component in place.
 *
 * Preserves field metadata and applies no finite-value validation.
 *
 * @param scalar Value added to each component.
 * @return Modified left operand.
 */
Field& Field::operator+=(Precision scalar) {
    Precision* field_data = values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        field_data[i] += scalar;
    }

    return *this;
}

/**
 * @brief Subtracts a scalar from every component in place.
 *
 * Preserves field metadata and applies no finite-value validation.
 *
 * @param scalar Value subtracted from each component.
 * @return Modified left operand.
 */
Field& Field::operator-=(Precision scalar) {
    Precision* field_data = values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        field_data[i] -= scalar;
    }

    return *this;
}

/**
 * @brief Multiplies every component by a scalar in place.
 *
 * This is a scalar-value operation, not a tensor product.
 *
 * @param scalar Uniform multiplier.
 * @return Modified left operand.
 */
Field& Field::operator*=(Precision scalar) {
    Precision* field_data = values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        field_data[i] *= scalar;
    }

    return *this;
}

/**
 * @brief Divides every component by a scalar in place.
 *
 * An exact zero divisor is rejected before any values are modified.
 * Nonfinite divisors and results are not otherwise checked.
 *
 * @param scalar Uniform divisor.
 * @return Modified left operand.
 */
Field& Field::operator/=(Precision scalar) {
    // Validate the uniform divisor before touching the destination.
    logging::error(scalar != Precision(0),
        "Field '", name, "': division by zero in scalar '/=' operation");

    Precision* field_data = values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        field_data[i] /= scalar;
    }

    return *this;
}

/**
 * @brief Adds a compatible field componentwise in place.
 *
 * Domain and dimensions must match; names, units and coordinate bases are
 * not compared. Self-use is supported through in-place arithmetic.
 *
 * @param other Field supplying the added components.
 * @return Modified left operand.
 */
Field& Field::operator+=(const Field& other) {
    // Confirm identical layouts before traversing the two scalar buffers.
    validate_compatible(other, "+=");

    Precision*       lhs = values.data();
    const Precision* rhs = other.values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        lhs[i] += rhs[i];
    }

    return *this;
}

/**
 * @brief Subtracts a compatible field componentwise in place.
 *
 * Domain and dimensions must match; names, units and coordinate bases are
 * not compared. Self-use is supported through in-place arithmetic.
 *
 * @param other Field supplying the subtracted components.
 * @return Modified left operand.
 */
Field& Field::operator-=(const Field& other) {
    // Confirm identical layouts before traversing the two scalar buffers.
    validate_compatible(other, "-=");

    Precision*       lhs = values.data();
    const Precision* rhs = other.values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        lhs[i] -= rhs[i];
    }

    return *this;
}

/**
 * @brief Multiplies compatible fields componentwise in place.
 *
 * This is an entrywise product, not matrix or tensor multiplication.
 * Domain and dimensions must match; physical meaning is caller-defined.
 *
 * @param other Field supplying the component multipliers.
 * @return Modified left operand.
 */
Field& Field::operator*=(const Field& other) {
    // Confirm identical layouts before traversing the two scalar buffers.
    validate_compatible(other, "*=");

    Precision*       lhs = values.data();
    const Precision* rhs = other.values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        lhs[i] *= rhs[i];
    }

    return *this;
}

/**
 * @brief Divides compatible fields componentwise in place.
 *
 * Domain and dimensions must match. Each denominator is checked for exact
 * zero immediately before division; a failure may leave earlier entries
 * modified. Units and component bases are not verified.
 *
 * @param other Field supplying the component divisors.
 * @return Modified left operand.
 */
Field& Field::operator/=(const Field& other) {
    // Validate layout once; denominators are checked during each update.
    validate_compatible(other, "/=");

    // Check each component so diagnostics identify its precise location.
    for (Index row = 0; row < rows; ++row) {
        for (Index component = 0; component < components; ++component) {
            // A zero in a later component leaves earlier divisions applied.
            const Precision denominator = other(row, component);
            logging::error(denominator != Precision(0),
                "Field '", name, "': division by zero at (", row, ",", component,
                ") in field '/=' with '", other.name, "'");

            values(row, component) /= denominator;
        }
    }

    return *this;
}

/**
 * @brief Returns the componentwise sum of two fields.
 *
 * Copies the left operand, then applies `+=`. The new field retains the left
 * operand's metadata, while both input fields remain unchanged.
 *
 * @param other Field added to this field.
 * @return Independent sum field.
 */
Field Field::operator+(const Field& other) const {
    // Compound arithmetic validates the inputs and preserves the left metadata.
    Field result = *this;
    result += other;
    return result;
}

/**
 * @brief Returns the componentwise difference of two fields.
 *
 * Copies the left operand, then applies `-=`. The new field retains the left
 * operand's metadata, while both input fields remain unchanged.
 *
 * @param other Field subtracted from this field.
 * @return Independent difference field.
 */
Field Field::operator-(const Field& other) const {
    // Compound arithmetic validates the inputs and preserves the left metadata.
    Field result = *this;
    result -= other;
    return result;
}

/**
 * @brief Returns the componentwise product of two fields.
 *
 * Copies the left operand, then applies `*=`. The operation is entrywise and
 * does not represent tensor or matrix multiplication.
 *
 * @param other Field multiplied with this field.
 * @return Independent product field.
 */
Field Field::operator*(const Field& other) const {
    // Compound arithmetic validates the inputs and preserves the left metadata.
    Field result = *this;
    result *= other;
    return result;
}

/**
 * @brief Returns the componentwise quotient of two fields.
 *
 * Copies the left operand, then applies `/=`. Denominator checks are inherited
 * from the compound operator; the input fields are not modified.
 *
 * @param other Field dividing this field.
 * @return Independent quotient field.
 */
Field Field::operator/(const Field& other) const {
    // Compound arithmetic validates the inputs and preserves the left metadata.
    Field result = *this;
    result /= other;
    return result;
}

/**
 * @brief Validates layouts before componentwise field operations.
 *
 * Requires identical domains, row counts and component counts. Names, units,
 * coordinate bases and physical component meanings are not compared.
 *
 * @param other Field whose metadata is compared.
 * @param operation Operator name used in an error message.
 */
void Field::validate_compatible(const Field& other,
                                const char*  operation) const {
    // Compare the entity domain separately from the rectangular field shape.
    logging::error(domain == other.domain,
        "Field '", name, "': domain mismatch in '", operation, "' (",
        static_cast<int>(domain), " vs ", static_cast<int>(other.domain), ")"
    );

    logging::error(rows == other.rows && components == other.components,
        "Field '", name, "': size mismatch in '", operation, "' (",
        rows, "x", components, " vs ", other.rows, "x", other.components, ")"
    );
}

} // namespace model
} // namespace fem

