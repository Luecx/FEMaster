#pragma once

#include "../core/types.h"
#include "csqrt.h"
#include <functional>
#include <ostream>
#include <map>
#include <vector>
#include <tuple>

namespace fem {

namespace quadrature {

enum Order {
    ORDER_CONSTANT        = 0,
    ORDER_LINEAR          = 1,
    ORDER_QUADRATIC       = 2,
    ORDER_CUBIC           = 3,
    ORDER_QUARTIC         = 4,
    ORDER_QUINTIC         = 5,
    ORDER_SUPER_LINEAR    = 1,
    ORDER_SUPER_QUADRATIC = 2,
};

enum Domain {
    DOMAIN_ISO_TRI,
    DOMAIN_ISO_QUAD,
    DOMAIN_ISO_HEX,
    DOMAIN_ISO_TET,
    DOMAIN_ISO_WEDGE,
    DOMAIN_ISO_LINE_A,
    DOMAIN_ISO_LINE_B,
    DOMAIN_ISO_PYRAMID,
};

struct Point {
    Precision r;
    Precision s;
    Precision t;
    Precision w;

    // Default constructor with member initialization
    constexpr Point() : r(0), s(0), t(0), w(0) {}

    // Constructor with all parameters
    constexpr Point(Precision r, Precision s, Precision t, Precision w)
        : r(r), s(s), t(t), w(w) {}

    // Constructor with three parameters (r, s, w), set t to 0
    constexpr Point(Precision r, Precision s, Precision w)
        : r(r), s(s), t(0), w(w) {}

    // Constructor with two parameters (r, w), set s and t to 0
    constexpr Point(Precision r, Precision w)
        : r(r), s(0), t(0), w(w) {}

    // Function to apply the weights
    template<typename T>
    constexpr T operator()(const std::function<T(Precision, Precision, Precision)>& func) const {
        return func(r, s, t) * w;
    }

    // Overloaded operator for printing
    friend std::ostream& operator<<(std::ostream& os, const Point& point) {
        os << "r: " << point.r << " s: " << point.s << " t: " << point.t << " w: " << point.w;
        return os;
    }
};

struct Scheme {
    Domain             domain;
    Order              order;
    std::vector<Point> points;

    // Default constructor with member initialization
    Scheme() : domain(Domain::DOMAIN_ISO_TRI), order(Order::ORDER_CONSTANT), points() {}
    Scheme(Domain domain, Order order, std::vector<Point> points)
        : domain(domain), order(order), points(points) {}
};

template<std::size_t N>
std::vector<Point> create_vector_from_array(const std::array<Point, N>& arr) {
    return std::vector<Point>(arr.begin(), arr.end());
}

// one for c-arrays
template<std::size_t N>
std::vector<Point> create_vector_from_array(const Point (&arr)[N]) {
    return std::vector<Point>(arr, arr + N);
}

// Type alias for the map key
using SchemeKey = std::tuple<Domain, Order>;

extern std::map<SchemeKey, Scheme> quadrature_scheme_map;

struct Executor {
    Executor(std::function<void()> f) {
        f();
    }
};

#define REGISTER_SCHEME_IN_PLACE(domain_, order_, ...) \
    static Executor executor_##domain_##order_([]() { \
        quadrature_scheme_map[{domain_, order_}] = Scheme(domain_, order_,std::vector<Point>{__VA_ARGS__}); \
    })
#define REGISTER_SCHEME(domain_, order_, points_) \
    static Executor executor_##domain_##order_([]() { \
        quadrature_scheme_map[{domain_, order_}] = Scheme(domain_, order_, create_vector_from_array(points_)); \
    })

struct Quadrature {
    private:
    Domain             domain;
    Order              order;
    Scheme             scheme;

    public:
    Quadrature(Domain d, Order o);

    template<typename F, typename T = std::invoke_result_t<F, Precision, Precision, Precision>>
    T integrate(F const& func) const {
        T res = T();
        if constexpr (!std::is_arithmetic_v<T>) {
            res.setZero();
        }
        for (Index i = 0; i < count(); i++) {
            res += get_point(i)(func);
        }
        return res;
    }

    Point get_point(ID n) const;
    ID    count() const;
};

}    // namespace quadrature

}    // namespace fem
