/**
 * @file quadrature.cpp
 * @brief Implements quadrature scheme lookup helpers.
 *
 * This translation unit owns the global scheme registry used by the
enumeration-based lookup API and provides the light-weight `Quadrature`
constructor/queries.
 *
 * @see src/math/quadrature.h
 */

#include "quadrature.h"

#include <stdexcept>

namespace fem {
namespace quadrature {

std::map<SchemeKey, Scheme> quadrature_scheme_map{};

Quadrature::Quadrature(Domain d, Order o)
    : domain(d), order(o) {
    auto it = quadrature_scheme_map.find({domain, order});
    if (it == quadrature_scheme_map.end()) {
        throw std::runtime_error("Quadrature scheme not found.");
    }
    scheme = it->second;
}

Point Quadrature::get_point(ID n) const {
    return scheme.points[static_cast<std::size_t>(n)];
}

Index Quadrature::count() const {
    return static_cast<Index>(scheme.points.size());
}

} // namespace quadrature
} // namespace fem

