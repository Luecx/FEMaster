/******************************************************************************
 * @file quadrature.h
 * @brief Declares Gauss-type quadrature schemes used throughout the solver.
 *
 * The API exposes compile-time registered integration schemes for various
 * isoparametric domains (lines, triangles, quads, tetrahedra, etc.) and
 * supplies convenience routines to execute numerical integration.
 *
 * @see src/math/quadrature.cpp
 * @see src/math/quadrature_iso_tri.cpp
 * @see src/math/quadrature_iso_hex.cpp
 ******************************************************************************/

#pragma once

#include "../core/assert.h"
#include "../core/types_eig.h"
#include "csqrt.h"

#include <array>
#include <functional>
#include <map>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <vector>

namespace fem {
namespace quadrature {

/******************************************************************************
 * @enum Order
 * @brief Enumerates polynomial orders supported by the quadrature collection.
 ******************************************************************************/
enum Order {
    ORDER_CONSTANT = 0,
    ORDER_LINEAR = 1,
    ORDER_QUADRATIC = 2,
    ORDER_CUBIC = 3,
    ORDER_QUARTIC = 4,
    ORDER_QUINTIC = 5,
    ORDER_SUPER_LINEAR = 1,
    ORDER_SUPER_QUADRATIC = 2,
};

/******************************************************************************
 * @enum Domain
 * @brief Enumerates reference domains covered by the integration schemes.
 ******************************************************************************/
enum Domain {
    DOMAIN_ISO_TRI = 1,
    DOMAIN_ISO_QUAD = 2,
    DOMAIN_ISO_HEX = 3,
    DOMAIN_ISO_TET = 4,
    DOMAIN_ISO_WEDGE = 5,
    DOMAIN_ISO_LINE_A = 6,
    DOMAIN_ISO_LINE_B = 7,
    DOMAIN_ISO_PYRAMID = 8,
};

/******************************************************************************
 * @struct Point
 * @brief Represents a quadrature point with coordinates and weight.
 ******************************************************************************/
struct Point {
    Precision r{0}; ///< First isoparametric coordinate.
    Precision s{0}; ///< Second isoparametric coordinate.
    Precision t{0}; ///< Third isoparametric coordinate.
    Precision w{0}; ///< Quadrature weight.

    /// Default constructor.
    constexpr Point() = default;

    /// Constructs a point with fully specified coordinates and weight.
    constexpr Point(Precision r_in, Precision s_in, Precision t_in, Precision w_in)
        : r(r_in), s(s_in), t(t_in), w(w_in) {}

    /// Constructs a point for 2D schemes (implicitly `t = 0`).
    constexpr Point(Precision r_in, Precision s_in, Precision w_in)
        : r(r_in), s(s_in), t(0), w(w_in) {}

    /// Constructs a point for 1D schemes (implicitly `s = t = 0`).
    constexpr Point(Precision r_in, Precision w_in)
        : r(r_in), s(0), t(0), w(w_in) {}

    /// Applies the quadrature weight to the provided function value.
    template<typename T>
    constexpr T operator()(const std::function<T(Precision, Precision, Precision)>& func) const {
        return func(r, s, t) * w;
    }

    /// Streams the coordinates and weight to an output stream.
    friend std::ostream& operator<<(std::ostream& os, const Point& point) {
        os << "r: " << point.r << " s: " << point.s << " t: " << point.t << " w: " << point.w;
        return os;
    }
};

/******************************************************************************
 * @struct Scheme
 * @brief Bundles all points belonging to a domain/order combination.
 ******************************************************************************/
struct Scheme {
    Domain domain{DOMAIN_ISO_TRI}; ///< Reference domain of the scheme.
    Order order{ORDER_CONSTANT};   ///< Supported polynomial order.
    std::vector<Point> points{};   ///< Quadrature points and weights.

    Scheme() = default;
    Scheme(Domain domain_in, Order order_in, std::vector<Point> points_in)
        : domain(domain_in), order(order_in), points(std::move(points_in)) {}
};

/// Alias for the key type that indexes all schemes.
using SchemeKey = std::tuple<Domain, Order>;

/******************************************************************************
 * @brief Converts an `std::array` of points into a `std::vector`.
 ******************************************************************************/
template<std::size_t N>
std::vector<Point> create_vector_from_array(const std::array<Point, N>& arr) {
    return std::vector<Point>(arr.begin(), arr.end());
}

/******************************************************************************
 * @brief Converts a C-style array of points into a `std::vector`.
 ******************************************************************************/
template<std::size_t N>
std::vector<Point> create_vector_from_array(const Point (&arr)[N]) {
    return std::vector<Point>(arr, arr + N);
}

/// Global lookup table populated by the registration macros.
extern std::map<SchemeKey, Scheme> quadrature_scheme_map;

/******************************************************************************
 * @brief Helper that executes a callable at static initialisation time.
 ******************************************************************************/
struct Executor {
    explicit Executor(std::function<void()> f) {
        f();
    }
};

/// Registers a scheme by directly passing the point list to the constructor.
#define REGISTER_SCHEME_IN_PLACE(domain_, order_, ...) \
    static ::fem::quadrature::Executor executor_##domain_##order_([]() { \
        ::fem::quadrature::quadrature_scheme_map[{domain_, order_}] = \
        ::fem::quadrature::Scheme(domain_, order_, std::vector<::fem::quadrature::Point>{__VA_ARGS__}); \
    })

/// Registers a scheme by forwarding an existing array through `create_vector_from_array`.
#define REGISTER_SCHEME(domain_, order_, points_) \
    static ::fem::quadrature::Executor executor_##domain_##order_([]() { \
        ::fem::quadrature::quadrature_scheme_map[{domain_, order_}] = \
        ::fem::quadrature::Scheme(domain_, order_, ::fem::quadrature::create_vector_from_array(points_)); \
    })

/******************************************************************************
 * @class Quadrature
 * @brief Provides access to a registered quadrature scheme.
 ******************************************************************************/
struct Quadrature {
private:
    Domain domain; ///< Domain associated with the instance.
    Order order;   ///< Order associated with the instance.
    Scheme scheme; ///< Cached scheme retrieved from the registry.

public:
    /******************************************************************************
     * @brief Constructs a quadrature instance for the specified domain/order pair.
     ******************************************************************************/
    Quadrature(Domain d, Order o);

    /******************************************************************************
     * @brief Integrates a callable using all points of the scheme.
     *
     * @tparam F Integrand functor type.
     * @tparam T Result type inferred from the callable.
     ******************************************************************************/
    template<typename F, typename T = std::invoke_result_t<F, Precision, Precision, Precision>>
    T integrate(const F& func) const {
        T res = T();
        if constexpr (!std::is_arithmetic_v<T>) {
            res.setZero();
        }
        for (Index i = 0; i < count(); ++i) {
            res += get_point(i)(func);
        }
        return res;
    }

    /******************************************************************************
     * @brief Extrapolates point-wise values back to nodal values (currently unimplemented).
     ******************************************************************************/
    template<int M, int N>
    StaticVector<N> extrapolate(const StaticMatrix<M, N>& values) const {
        runtime_assert(M == count(), "Number of values must match number of quadrature points");
        return StaticVector<N>::Zero();
    }

    /// Returns the `n`-th quadrature point.
    Point get_point(ID n) const;

    /// Returns the number of quadrature points in the scheme.
    Index count() const;
};

} // namespace quadrature
} // namespace fem

