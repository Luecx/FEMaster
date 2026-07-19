/**
 * @file field.cpp
 * @brief Implements lightweight dense storage for model fields.
 *
 * @see field.h
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

FieldMatrix::FieldMatrix(Index rows, Index cols)
    : rows_(rows),
      cols_(cols) {
    logging::error(rows >= 0, "FieldMatrix: rows must be non-negative");
    logging::error(cols >= 0, "FieldMatrix: cols must be non-negative");

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

bool FieldMatrix::has_any_finite() const {
    for (const Precision value : data_) {
        if (std::isfinite(static_cast<double>(value))) {
            return true;
        }
    }

    return false;
}

std::size_t FieldMatrix::offset(Index row, Index col) const {
    logging::error(row >= 0 && row < rows_ && col >= 0 && col < cols_,
        "FieldMatrix: index (", row, ", ", col, ") is outside ", rows_, "x", cols_);

    return static_cast<std::size_t>(row) *
           static_cast<std::size_t>(cols_) +
           static_cast<std::size_t>(col);
}

// ------------------------------------------------------------
// Field
// ------------------------------------------------------------

Field::Field(std::string field_name,
             FieldDomain field_domain,
             Index       row_count,
             Index       component_count)
    : name      (std::move(field_name)),
      domain    (field_domain),
      rows      (row_count),
      components(component_count),
      values    (rows, components) {
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

Precision& Field::operator()(Index row) {
    logging::error(components == 1,
        "Field '", name, "': scalar access requires exactly one component");

    return values(row, 0);
}

Precision Field::operator()(Index row) const {
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

bool Field::is_nan(Index row, Index component) const {
    const Precision value = values(row, component);

    return !std::isfinite(static_cast<double>(value));
}

Vec3 Field::row_vec3(Index row) const {
    logging::error(components >= 3,
        "Field '", name, "': row_vec3 requires at least three components");

    return Vec3(
        values(row, 0),
        values(row, 1),
        values(row, 2)
    );
}

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

Field& Field::operator+=(Precision scalar) {
    Precision* field_data = values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        field_data[i] += scalar;
    }

    return *this;
}

Field& Field::operator-=(Precision scalar) {
    Precision* field_data = values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        field_data[i] -= scalar;
    }

    return *this;
}

Field& Field::operator*=(Precision scalar) {
    Precision* field_data = values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        field_data[i] *= scalar;
    }

    return *this;
}

Field& Field::operator/=(Precision scalar) {
    logging::error(scalar != Precision(0),
        "Field '", name, "': division by zero in scalar '/=' operation");

    Precision* field_data = values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        field_data[i] /= scalar;
    }

    return *this;
}

Field& Field::operator+=(const Field& other) {
    validate_compatible(other, "+=");

    Precision*       lhs = values.data();
    const Precision* rhs = other.values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        lhs[i] += rhs[i];
    }

    return *this;
}

Field& Field::operator-=(const Field& other) {
    validate_compatible(other, "-=");

    Precision*       lhs = values.data();
    const Precision* rhs = other.values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        lhs[i] -= rhs[i];
    }

    return *this;
}

Field& Field::operator*=(const Field& other) {
    validate_compatible(other, "*=");

    Precision*       lhs = values.data();
    const Precision* rhs = other.values.data();

    for (std::size_t i = 0; i < values.size(); ++i) {
        lhs[i] *= rhs[i];
    }

    return *this;
}

Field& Field::operator/=(const Field& other) {
    validate_compatible(other, "/=");

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

void Field::validate_compatible(const Field& other,
                                const char*  operation) const {
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
