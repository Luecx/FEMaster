#pragma once

#include "../core/types.h"

#include <functional>
#include <ostream>

namespace fem {

namespace quadrature {

enum Order {
    ORDER_CONSTANT        = 0,
    ORDER_LINEAR          = 1,
    ORDER_QUADRATIC       = 2,
    ORDER_CUBIC           = 3,
    ORDER_QUARTIC         = 4,
    ORDER_QUINTIC         = 5,
    ORDER_SUPER_LINEAR    = 1,    // used for wedge, linear on triangular side, quadratic along edge
    ORDER_SUPER_QUADRATIC = 2,    // used for wedge, quadratic on triangular side, quartic along edge
};

enum Domain {
    DOMAIN_ISO_TRI,
    DOMAIN_ISO_QUAD,
    DOMAIN_ISO_HEX,
    DOMAIN_ISO_TET,
    DOMAIN_ISO_WEDGE,
};

struct Point {
    Precision r;
    Precision s;
    Precision t;

    Precision w;
    Point() {}
    Point(Precision r, Precision s, Precision t, Precision w);
    Point(Precision r, Precision s, Precision w);

    template<typename T>
    T operator()(const std::function<T(Precision, Precision, Precision)>& func) const {
        return func(r, s, t) * w;
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& point) {
        os << "r: " << point.r << " s: " << point.s << " t: " << point.t << " w: " << point.w;
        return os;
    }
};

struct Quadrature {
    private:
    Domain             domain;
    Order              order;
    std::vector<Point> points {};

    public:
    Quadrature(Domain d, Order o);

    template<typename T>
    T integrate(const std::function<T(Precision, Precision, Precision)>& func) const {
        T res = T();
        if constexpr (!std::is_arithmetic_v<T>) {
            res.setZero();
        }
        for (int i = 0; i < points.size(); i++) {
            res += points[i](func);
        }
        return res;
    }

    Point get_point(ID n) const;

    ID    count() const;
};
}    // namespace quadrature

}    // namespace fem
