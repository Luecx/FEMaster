/**
 * @file field.cpp
 * @brief Implements dense field storage, validation and componentwise arithmetic.
 *
 * The model-data implementation allocates row-major scalar vectors, converts
 * checked row/component indices into storage offsets and exposes non-owning data
 * access. Field operations preserve metadata and update values in place after
 * checking the domain and dimensions where a second field participates.
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

// ------------------------------------------------------------
// FieldMatrix
// ------------------------------------------------------------

/**
 * @brief Allocates zero-initialized row-major matrix storage.
 *
 * Convert the unsigned dimensions to storage-size values and allocate rows * cols
 * Precision entries. Zero dimensions are permitted. The product must be
 * representable in std::size_t and fit available memory; overflow is not checked.
 *
 * @param rows Number of stored rows.
 * @param cols Number of scalar columns in each row.
 */
FieldMatrix::FieldMatrix(Index rows, Index cols)
    : rows_(rows),
      cols_(cols) {
    // Convert dimension metadata to vector sizes before allocating zero-initialized values.
    const auto row_count = static_cast<std::size_t>(rows);
    const auto col_count = static_cast<std::size_t>(cols);

    data_.resize(row_count * col_count);
}

Index FieldMatrix::rows() const {
    return rows_;
}

Index FieldMatrix::cols() const {
    return cols_;
}

std::size_t FieldMatrix::size() const {
    return data_.size();
}

Precision& FieldMatrix::operator()(Index row, Index col) {
    return data_[offset(row, col)];
}

Precision FieldMatrix::operator()(Index row, Index col) const {
    return data_[offset(row, col)];
}

Precision* FieldMatrix::data() {
    return data_.data();
}

const Precision* FieldMatrix::data() const {
    return data_.data();
}

void FieldMatrix::set_zero() {
    std::fill(data_.begin(), data_.end(), Precision(0));
}

void FieldMatrix::set_ones() {
    std::fill(data_.begin(), data_.end(), Precision(1));
}

void FieldMatrix::fill_nan() {
    const Precision nan = std::numeric_limits<Precision>::quiet_NaN();

    std::fill(data_.begin(), data_.end(), nan);
}

/**
 * @brief Tests whether dense storage contains at least one finite scalar.
 *
 * Scan values until std::isfinite succeeds, converting Precision to double for
 * the predicate. Empty or entirely NaN/infinite storage returns false. Values and
 * dimensions remain unchanged; this query does not require all entries to be valid.
 *
 * @return True if at least one stored entry is finite.
 */
bool FieldMatrix::has_any_finite() const {
    // Stop at the first usable scalar; empty and all-nonfinite matrices have none.
    for (const Precision value : data_) {
        if (std::isfinite(static_cast<double>(value))) {
            return true;
        }
    }

    return false;
}

/**
 * @brief Checks matrix bounds and computes the row-major storage offset.
 *
 * Require row < rows_ and col < cols_ before evaluating row * cols_ + col.
 * Indices are unsigned, so negative input is not a valid caller convention.
 * Storage dimensions must have a representable product. The check reads metadata
 * and does not change matrix storage.
 *
 * @param row Zero-based row index.
 * @param col Zero-based scalar column index.
 * @return Flat index into the owned row-major value vector.
 */
std::size_t FieldMatrix::offset(Index row, Index col) const {
    // Check both component bounds before converting the matrix index into a flat offset.
    logging::error(row < rows_ && col < cols_,
        "FieldMatrix: index (", row, ", ", col, ") is outside ", rows_, "x", cols_);

    return static_cast<std::size_t>(row) *
           static_cast<std::size_t>(cols_) +
           static_cast<std::size_t>(col);
}

// ------------------------------------------------------------
// Field
// ------------------------------------------------------------

/**
 * @brief Constructs named domain metadata and zero-initialized dense values.
 *
 * Initialize public metadata and allocate FieldMatrix(rows, components), then
 * check that both dimensions are positive. Allocation precedes these checks;
 * callers must supply a representable storage-size product. No row remapping,
 * physical unit assignment or basis transformation is performed.
 *
 * @param field_name Mutable identifier used in field diagnostics.
 * @param field_domain Model entity layout associated with each row.
 * @param row_count Positive number of rows in the compiled domain layout.
 * @param component_count Positive number of scalar components per row.
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
    // Validate positive public dimensions after the owned values have been initialized.
    logging::error(rows > 0,
        "Field '", name, "': rows must be positive");
    logging::error(components > 0,
        "Field '", name, "': components must be positive");
}

Precision& Field::operator()(Index row, Index component) {
    return values(row, component);
}

Precision Field::operator()(Index row, Index component) const {
    return values(row, component);
}

/**
 * @brief Accesses the sole component of a scalar field at a checked row.
 *
 * Require exactly one component before delegating to matrix indexing. This form
 * selects scalar storage without changing metadata or interpreting physical units.
 *
 * @param row Zero-based field row.
 * @return Mutable reference to the scalar entry.
 */
Precision& Field::operator()(Index row) {
    // Reject ambiguous scalar access before selecting column zero.
    logging::error(components == 1,
        "Field '", name, "': scalar access requires exactly one component");

    return values(row, 0);
}

/**
 * @brief Accesses the sole component of a scalar field at a checked row.
 *
 * Require exactly one component before delegating to matrix indexing. This form
 * selects scalar storage without changing metadata or interpreting physical units.
 *
 * @param row Zero-based field row.
 * @return Copy of the scalar entry.
 */
Precision Field::operator()(Index row) const {
    // Reject ambiguous scalar access before selecting column zero.
    logging::error(components == 1,
        "Field '", name, "': scalar access requires exactly one component");

    return values(row, 0);
}

Precision* Field::data() {
    return values.data();
}

const Precision* Field::data() const {
    return values.data();
}

void Field::set_zero() {
    values.set_zero();
}

void Field::set_ones() {
    values.set_ones();
}

void Field::fill_nan() {
    values.fill_nan();
}

bool Field::has_any_finite() const {
    return values.has_any_finite();
}

/**
 * @brief Checks a stored component for any nonfinite value.
 *
 * Read the entry through bounds-checked matrix access and negate std::isfinite.
 * Both NaN and positive/negative infinity satisfy this predicate; the operation
 * does not modify field storage or metadata.
 *
 * @param row Zero-based field row.
 * @param component Zero-based scalar component.
 * @return True for NaN or infinity, false for a finite scalar.
 */
bool Field::is_nan(Index row, Index component) const {
    // Use the matrix bounds check before testing the complete nonfinite category.
    const Precision value = values(row, component);

    return !std::isfinite(static_cast<double>(value));
}

/**
 * @brief Validates that every stored field component is finite.
 *
 * NaN and infinite values are treated as numerical failures. The supplied
 * label identifies the operation or physical result being checked, while the
 * reported row and component locate the invalid value in the field storage.
 * The method reads rows/components metadata, which must match values storage,
 * and leaves the field unchanged. An empty default field has no entries to check.
 *
 * @param label Diagnostic name included in a failed check.
 */
void Field::check_finite(const std::string& label) const {
    // Inspect the complete rectangular field storage
    for (Index row = 0; row < rows; ++row) {
        for (Index component = 0; component < components; ++component) {
            logging::error(std::isfinite(static_cast<double>(values(row, component))),
                label, " row ", row, " has invalid value at component ", component);
        }
    }
}

/**
 * @brief Copies the first 3 components of one field row into Vec3.
 *
 * Require at least 3 components and read them in stored order through checked
 * matrix access. Additional components are ignored. No tensor-component
 * reordering, physical interpretation or coordinate transformation is performed.
 *
 * @param row Zero-based field row.
 * @return Independent vector copy in the caller-defined component basis.
 */
Vec3 Field::row_vec3(Index row) const {
    // Ensure the row has enough components before copying its first 3 scalar values.
    logging::error(components >= 3,
        "Field '", name, "': row_vec3 requires at least three components");

    return Vec3(
        values(row, 0),
        values(row, 1),
        values(row, 2)
    );
}

/**
 * @brief Copies the first 6 components of one field row into Vec6.
 *
 * Require at least 6 components and read them in stored order through checked
 * matrix access. Additional components are ignored. No tensor-component
 * reordering, physical interpretation or coordinate transformation is performed.
 *
 * @param row Zero-based field row.
 * @return Independent vector copy in the caller-defined component basis.
 */
Vec6 Field::row_vec6(Index row) const {
    // Ensure the row has enough components before copying its first 6 scalar values.
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
 * @brief Applies scalar addition to every stored field component.
 *
 * Update the contiguous row-major values with a_i <- a_i + scalar, retaining the name,
 * domain and dimensions. This is componentwise arithmetic in the existing
 * caller-defined basis; there is no vector/tensor transformation.
 * Nonfinite inputs/results are not checked by this operation.
 *
 * @param scalar Scalar operand applied uniformly to all entries.
 * @return This field after the in-place update.
 */
Field& Field::operator+=(Precision scalar) {
    // Traverse owned contiguous storage without changing the field metadata.
    Precision* field_data = values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        field_data[i] += scalar;
    }

    return *this;
}

/**
 * @brief Applies scalar subtraction to every stored field component.
 *
 * Update the contiguous row-major values with a_i <- a_i - scalar, retaining the name,
 * domain and dimensions. This is componentwise arithmetic in the existing
 * caller-defined basis; there is no vector/tensor transformation.
 * Nonfinite inputs/results are not checked by this operation.
 *
 * @param scalar Scalar operand applied uniformly to all entries.
 * @return This field after the in-place update.
 */
Field& Field::operator-=(Precision scalar) {
    // Traverse owned contiguous storage without changing the field metadata.
    Precision* field_data = values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        field_data[i] -= scalar;
    }

    return *this;
}

/**
 * @brief Applies scalar multiplication to every stored field component.
 *
 * Update the contiguous row-major values with a_i <- a_i * scalar, retaining the name,
 * domain and dimensions. This is componentwise arithmetic in the existing
 * caller-defined basis; there is no vector/tensor transformation.
 * Nonfinite inputs/results are not checked by this operation.
 *
 * @param scalar Scalar operand applied uniformly to all entries.
 * @return This field after the in-place update.
 */
Field& Field::operator*=(Precision scalar) {
    // Traverse owned contiguous storage without changing the field metadata.
    Precision* field_data = values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        field_data[i] *= scalar;
    }

    return *this;
}

/**
 * @brief Applies scalar division to every stored field component.
 *
 * Update the contiguous row-major values with a_i <- a_i / scalar, retaining the name,
 * domain and dimensions. This is componentwise arithmetic in the existing
 * caller-defined basis; there is no vector/tensor transformation.
 * Require scalar != 0 before any update; nonfinite inputs/results are not checked.
 *
 * @param scalar Scalar operand applied uniformly to all entries.
 * @return This field after the in-place update.
 */
Field& Field::operator/=(Precision scalar) {
    // Reject an exact zero scalar before modifying any component.
    logging::error(scalar != Precision(0),
        "Field '", name, "': division by zero in scalar '/=' operation");

    // Traverse owned contiguous storage without changing the field metadata.
    Precision* field_data = values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        field_data[i] /= scalar;
    }

    return *this;
}

/**
 * @brief Performs componentwise field addition in the existing storage basis.
 *
 * Validate equal domains and dimensions before pairing entries with the same
 * row/component index. Names, physical units, tensor meaning and coordinate bases
 * are not compared; their compatibility is a caller responsibility. Update this
 * field in place and retain its metadata. The operand is not copied, so self-use
 * follows the same componentwise arithmetic.
 * Nonfinite values/results are not checked here.
 *
 * @param other Operand with matching row layout, domain and component count.
 * @return This field after the in-place component updates.
 */
Field& Field::operator+=(const Field& other) {
    // Establish matching row/component layouts before touching destination values.
    validate_compatible(other, "+=");

    // Pair entries by their shared row-major offset in the existing component basis.
    Precision*       lhs = values.data();
    const Precision* rhs = other.values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        lhs[i] += rhs[i];
    }

    return *this;
}

/**
 * @brief Performs componentwise field subtraction in the existing storage basis.
 *
 * Validate equal domains and dimensions before pairing entries with the same
 * row/component index. Names, physical units, tensor meaning and coordinate bases
 * are not compared; their compatibility is a caller responsibility. Update this
 * field in place and retain its metadata. The operand is not copied, so self-use
 * follows the same componentwise arithmetic.
 * Nonfinite values/results are not checked here.
 *
 * @param other Operand with matching row layout, domain and component count.
 * @return This field after the in-place component updates.
 */
Field& Field::operator-=(const Field& other) {
    // Establish matching row/component layouts before touching destination values.
    validate_compatible(other, "-=");

    // Pair entries by their shared row-major offset in the existing component basis.
    Precision*       lhs = values.data();
    const Precision* rhs = other.values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        lhs[i] -= rhs[i];
    }

    return *this;
}

/**
 * @brief Performs componentwise field multiplication in the existing storage basis.
 *
 * Validate equal domains and dimensions before pairing entries with the same
 * row/component index. Names, physical units, tensor meaning and coordinate bases
 * are not compared; their compatibility is a caller responsibility. Update this
 * field in place and retain its metadata. The operand is not copied, so self-use
 * follows the same componentwise arithmetic.
 * Multiplication, when used, is a scalar-entry product rather than a tensor product.
 * Nonfinite values/results are not checked here.
 *
 * @param other Operand with matching row layout, domain and component count.
 * @return This field after the in-place component updates.
 */
Field& Field::operator*=(const Field& other) {
    // Establish matching row/component layouts before touching destination values.
    validate_compatible(other, "*=");

    // Pair entries by their shared row-major offset in the existing component basis.
    Precision*       lhs = values.data();
    const Precision* rhs = other.values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        lhs[i] *= rhs[i];
    }

    return *this;
}

/**
 * @brief Performs componentwise field division in the existing storage basis.
 *
 * Validate equal domains and dimensions before pairing entries with the same
 * row/component index. Names, physical units, tensor meaning and coordinate bases
 * are not compared; their compatibility is a caller responsibility. Update this
 * field in place and retain its metadata. The operand is not copied, so self-use
 * follows the same componentwise arithmetic.
 * Each denominator is checked for exact zero immediately before its update.
 * A later failure can leave earlier components modified; there is no rollback.
 * Nonfinite values/results are not checked here.
 *
 * @param other Operand with matching row layout, domain and component count.
 * @return This field after the in-place component updates.
 */
Field& Field::operator/=(const Field& other) {
    // Establish matching row/component layouts before touching destination values.
    validate_compatible(other, "/=");

    // Check and divide one component at a time; no rollback follows a later zero divisor.
    for (Index row = 0; row < rows; ++row) {
        for (Index component = 0; component < components; ++component) {
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
 * @brief Validates structural compatibility for componentwise field arithmetic.
 *
 * Require the same FieldDomain, row count and component count before an operation
 * can pair stored entries. Public metadata must agree with each field's matrix
 * storage; this method does not check that invariant or physical units/bases.
 * No values or metadata are changed by validation.
 *
 * @param other Field whose structural metadata is compared with this field.
 * @param operation Diagnostic operator label included in failure messages.
 */
void Field::validate_compatible(const Field& other,
                                const char*  operation) const {
    // Verify semantic row association and rectangular shape as one compatibility phase.
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
