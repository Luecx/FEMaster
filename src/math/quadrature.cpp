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

SchemeRegistry& quadrature_scheme_registry() {
    static SchemeRegistry registry{};
    return registry;
}

void register_scheme(Domain domain, Order order, std::vector<Point> points) {
    quadrature_scheme_registry()[{domain, order}] = Scheme(domain, order, std::move(points));
}

Quadrature::Quadrature(Domain d, Order o)
    : domain(d), order(o) {
    auto& registry = quadrature_scheme_registry();
    auto it = registry.find({domain, order});
    if (it == registry.end()) {
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
