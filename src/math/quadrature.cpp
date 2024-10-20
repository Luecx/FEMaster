#include "quadrature.h"
#include "../core/assert.h"
#include <cassert>

namespace fem {
namespace quadrature {

// Initialize the global map
std::map<SchemeKey, Scheme> quadrature_scheme_map;

Quadrature::Quadrature(Domain d, Order o) : domain(d), order(o) {
    auto it = quadrature_scheme_map.find({domain, order});
    if (it != quadrature_scheme_map.end()) {
        scheme = it->second;
    } else {
        throw std::runtime_error("Quadrature scheme not found.");
    }
}

Point Quadrature::get_point(ID n) const {
    return scheme.points[n];
}

ID Quadrature::count() const {
    return scheme.points.size();
}

} // namespace quadrature
} // namespace fem
