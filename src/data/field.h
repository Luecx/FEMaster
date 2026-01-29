/**
 * @file field.h
 * @brief Declares a lightweight dense field container used for model data.
 *
 * Fields are intentionally untyped: they only store their name, domain, size,
 * and values. Semantics are handled by higher-level model utilities.
 *
 * @see src/model/model_data.h
 */

#pragma once

#include "../core/logging.h"
#include "../core/types_eig.h"

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace fem {
namespace model {

enum class FieldDomain : std::uint8_t {
    NODE,
    ELEMENT,
    IP
};

// ------------------------------------------------------------
// Dense storage (row-major): data[r*cols + c]
// ------------------------------------------------------------
class FieldMatrix {
public:
    FieldMatrix() = default;

    FieldMatrix(Index rows, Index cols)
        : _rows(rows), _cols(cols),
          _data(static_cast<size_t>(rows) * static_cast<size_t>(cols)) {}

    [[nodiscard]] Index rows() const { return _rows; }
    [[nodiscard]] Index cols() const { return _cols; }

    Precision& operator()(Index r, Index c) {
        return _data[static_cast<size_t>(r) * static_cast<size_t>(_cols) + static_cast<size_t>(c)];
    }

    Precision operator()(Index r, Index c) const {
        return _data[static_cast<size_t>(r) * static_cast<size_t>(_cols) + static_cast<size_t>(c)];
    }

    Precision* data() { return _data.data(); }
    const Precision* data() const { return _data.data(); }

    void set_zero() {
        for (auto& v : _data) v = Precision(0);
    }

    void set_ones() {
        for (auto& v : _data) v = Precision(1);
    }

    void fill_nan() {
        const Precision nan = std::numeric_limits<Precision>::quiet_NaN();
        for (auto& v : _data) v = nan;
    }

    [[nodiscard]] bool has_any_finite() const {
        for (const auto& v : _data) {
            if (std::isfinite(static_cast<double>(v))) return true;
        }
        return false;
    }

private:
    Index _rows = 0;
    Index _cols = 0;
    std::vector<Precision> _data{};
};

// ------------------------------------------------------------
// Field: name + domain + allocated sizes + dense data
// ------------------------------------------------------------
struct Field {
    using Ptr = std::shared_ptr<Field>;

    const std::string name;
    const FieldDomain domain;

    const Index rows;
    const Index components;

    FieldMatrix values;

    Field(std::string p_name, FieldDomain p_domain, Index p_rows, Index p_components)
        : name(std::move(p_name)),
          domain(p_domain),
          rows(p_rows),
          components(p_components),
          values(rows, components) {
        logging::error(rows > 0, "Field '", name, "': rows == 0 (domain not configured?)");
        logging::error(components > 0, "Field '", name, "': components == 0");
    }

    // ---- basic access ----
    Precision& operator()(Index r, Index c) { return values(r, c); }
    Precision  operator()(Index r, Index c) const { return values(r, c); }

    // convenience for scalar fields (components==1)
    Precision& operator()(Index r) {
        logging::error(components == 1, "Field '", name, "': scalar access requires components==1");
        return values(r, 0);
    }
    Precision operator()(Index r) const {
        logging::error(components == 1, "Field '", name, "': scalar access requires components==1");
        return values(r, 0);
    }

    Precision* data() { return values.data(); }
    const Precision* data() const { return values.data(); }

    // ---- init helpers ----
    void fill_nan() { values.fill_nan(); }
    void set_zero() { values.set_zero(); }
    void setOnes() { values.set_ones(); }

    // ---- validity helpers ----
    [[nodiscard]] bool has_any_finite() const { return values.has_any_finite(); }

    [[nodiscard]] bool is_nan(Index r, Index c) const {
        const Precision v = values(r, c);
        return !std::isfinite(static_cast<double>(v));
    }

    // ---- row helpers ----
    Vec3 row_vec3(Index r) const {
        logging::error(components >= 3, "Field '", name, "': row_vec3 requires components>=3");
        return Vec3(values(r, 0), values(r, 1), values(r, 2));
    }

    Vec6 row_vec6(Index r) const {
        logging::error(components >= 6, "Field '", name, "': row_vec6 requires components>=6");
        return Vec6(values(r, 0), values(r, 1), values(r, 2),
                    values(r, 3), values(r, 4), values(r, 5));
    }
};

} // namespace model
} // namespace fem
