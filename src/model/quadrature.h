#pragma once

#include "../core/types.h"

#include <functional>
#include <ostream>

namespace fem{

namespace quadrature{

enum Order {
    ORDER_LINEAR = 0,
    ORDER_QUADRATIC = 1,
    ORDER_CUBIC = 2,
    ORDER_BILINEAR = ORDER_QUADRATIC,
    ORDER_TRILINEAR = ORDER_CUBIC
};

enum Domain {
    DOMAIN_ISO_TRI = 3,
    DOMAIN_ISO_QUAD = 4,
    DOMAIN_ISO_HEX = 8
};

constexpr size_t n_dims(Order O, Domain D) {
    return D == DOMAIN_ISO_TRI ? 2 : D == DOMAIN_ISO_QUAD ? 2 : D == DOMAIN_ISO_HEX ? 3 : 3;
}

constexpr size_t n_points(Order O, Domain D) {
    return   (D == DOMAIN_ISO_TRI)  ?((O == ORDER_LINEAR)    ? 1
                                    : (O == ORDER_QUADRATIC) ? 3
                                    : (O == ORDER_CUBIC    ) ? 4
                                                             : 0)
           : (D == DOMAIN_ISO_QUAD) ?((O == ORDER_LINEAR)    ? 1
                                    : (O == ORDER_QUADRATIC) ? 4
                                    : (O == ORDER_CUBIC    ) ? 9
                                                             : 0)
           : (D == DOMAIN_ISO_HEX)  ?((O == ORDER_LINEAR)    ? 1
                                    : (O == ORDER_QUADRATIC) ? 8
                                    : (O == ORDER_CUBIC    ) ? 27
                                                             : 0)
           : 0;
}

struct Point {
    Precision r;
    Precision s;
    Precision t;

    Precision w;
    Point() {}
    Point(Precision r, Precision s, Precision t, Precision w)
        : r(r)
        , s(s)
        , t(t)
        , w(w) {}
    Point(Precision r, Precision s, Precision w)
        : r(r)
        , s(s)
        , t(0)
        , w(w) {}

    template<typename T>
    T operator()(const std::function<T(Precision, Precision, Precision)>& func) const{
        return func(r, s, t) * w;
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& point) {
        os << "r: " << point.r << " s: " << point.s << " t: " << point.t << " w: " << point.w;
        return os;
    }
};


struct Quadrature{
    private:
    Domain domain;
    Order order;
    std::vector<Point> points{};

    public:
    Quadrature(Domain d, Order o) : domain(d), order(o) {
        points.resize(n_points(o, d));

        if (d == DOMAIN_ISO_TRI){
            if (o == ORDER_LINEAR){
                points = {Point(1.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0)};
            }

            if (o == ORDER_QUADRATIC){
                points[0] = Point(1.0/6.0, 1.0/6.0, 1.0/6.0),
                points[1] = Point(2.0/3.0, 1.0/6.0, 1.0/6.0),
                points[2] = Point(1.0/6.0, 2.0/3.0, 1.0/6.0);
            }

            if (o == ORDER_CUBIC){
                points[0] = Point(1.0/3.0, 1.0/3.0, -27.0/48.0);
                points[1] = Point(1.0/5.0, 3.0/5.0,  25.0/48.0);
                points[2] = Point(1.0/5.0, 1.0/5.0,  25.0/48.0);
                points[4] = Point(3.0/5.0, 1.0/5.0,  25.0/48.0);
            }
        }

        if (d == DOMAIN_ISO_QUAD){

            if (o == ORDER_LINEAR){
                points[0] = Point(0.0, 0.0, 0.0, 4.0);
            }

            if (o == ORDER_QUADRATIC){
                points[0] = Point( sqrt(3.0)/3.0,  sqrt(3.0)/3.0, 1.0);
                points[1] = Point(-sqrt(3.0)/3.0,  sqrt(3.0)/3.0, 1.0);
                points[2] = Point(-sqrt(3.0)/3.0, -sqrt(3.0)/3.0, 1.0);
                points[3] = Point( sqrt(3.0)/3.0, -sqrt(3.0)/3.0, 1.0);
            }

            if (o == ORDER_CUBIC){
                points[0] = Point(-sqrt(0.6), -sqrt(0.6), 25.0 / 81.0);
                points[1] = Point( sqrt(0.0), -sqrt(0.6), 40.0 / 81.0);
                points[2] = Point( sqrt(0.6), -sqrt(0.6), 25.0 / 81.0);
                points[3] = Point(-sqrt(0.6),  sqrt(0.0), 40.0 / 81.0);
                points[4] = Point( sqrt(0.0),  sqrt(0.0), 64.0 / 81.0);
                points[5] = Point( sqrt(0.6),  sqrt(0.0), 40.0 / 81.0);
                points[6] = Point(-sqrt(0.6),  sqrt(0.6), 25.0 / 81.0);
                points[7] = Point( sqrt(0.0),  sqrt(0.6), 40.0 / 81.0);
                points[8] = Point( sqrt(0.6),  sqrt(0.6), 25.0 / 81.0);
            }
        }


        if (d == DOMAIN_ISO_HEX){
            if (o == ORDER_LINEAR){
                points[0] = Point(0.0, 0.0, 0.0, 8.0);
            }

            if (o == ORDER_QUADRATIC){
                points[0] = Point( 1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0), 1.0);
                points[1] = Point( 1.0/sqrt(3.0),  1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0);
                points[2] = Point( 1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0), 1.0);
                points[3] = Point( 1.0/sqrt(3.0), -1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0);
                points[4] = Point(-1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0), 1.0);
                points[5] = Point(-1.0/sqrt(3.0),  1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0);
                points[6] = Point(-1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0), 1.0);
                points[7] = Point(-1.0/sqrt(3.0), -1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0);
            }

            if (o == ORDER_CUBIC){
                Precision coords[3] {(Precision)-sqrt(3.0 / 5.0), 0, (Precision)sqrt(3.0 / 5.0)};
                Precision weight[3] { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
                int c = 0;
                for (int x = 0; x < 3; x++){
                    for (int y = 0; y < 3; y++){
                        for (int z = 0; z < 3; z++){
                            points[c] = Point(coords[x], coords[y], coords[z], weight[x] * weight[y] * weight[z]);
                            c ++;
                        }
                    }
                }
            }
        }

    }

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
};
}    // namespace quadrature

}    // namespace fem
