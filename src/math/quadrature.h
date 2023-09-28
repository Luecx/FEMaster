#pragma once

#include "../core/types.h"

#include <functional>
#include <ostream>

namespace fem{

namespace quadrature{

enum Order {
    ORDER_CONSTANT = 0,
    ORDER_LINEAR = 1,
    ORDER_QUADRATIC = 2,
    ORDER_CUBIC = 3,
    ORDER_QUARTIC = 4,
    ORDER_QUINTIC = 5,
    ORDER_SUPER_LINEAR = 1, // used for wedge, linear on triangular side, quadratic along edge
    ORDER_SUPER_QUADRATIC = 2, // used for wedge, quadratic on triangular side, quartic along edge
};

enum Domain {
    DOMAIN_ISO_TRI,
    DOMAIN_ISO_QUAD,
    DOMAIN_ISO_HEX,
    DOMAIN_ISO_TET,
    DOMAIN_ISO_WEDGE,
};

constexpr size_t n_dims(Order O, Domain D) {
    return D == DOMAIN_ISO_TRI ? 2
         : D == DOMAIN_ISO_QUAD ? 2
         : D == DOMAIN_ISO_HEX ? 3
         : D == DOMAIN_ISO_TET ? 3
         : D == DOMAIN_ISO_WEDGE ? 3
         : 3;

}

constexpr size_t n_points(Order O, Domain D) {
    return   (D == DOMAIN_ISO_TRI)  ?((O == ORDER_CONSTANT       ) ? 1
                                    : (O == ORDER_LINEAR         ) ? 1
                                    : (O == ORDER_QUADRATIC      ) ? 3
                                    : (O == ORDER_CUBIC          ) ? 4
                                                                   : 0)
           : (D == DOMAIN_ISO_QUAD) ?((O == ORDER_CONSTANT       ) ? 1
                                    : (O == ORDER_LINEAR         ) ? 1
                                    : (O == ORDER_QUADRATIC      ) ? 4
                                    : (O == ORDER_CUBIC          ) ? 4
                                    : (O == ORDER_QUARTIC        ) ? 9
                                    : (O == ORDER_QUINTIC        ) ? 9
                                                                   : 0)
           : (D == DOMAIN_ISO_HEX)  ?((O == ORDER_CONSTANT       ) ? 1
                                    : (O == ORDER_LINEAR         ) ? 1
                                    : (O == ORDER_QUADRATIC      ) ? 8
                                    : (O == ORDER_CUBIC          ) ? 8
                                    : (O == ORDER_QUARTIC        ) ? 27
                                    : (O == ORDER_QUINTIC        ) ? 27
                                                                   : 0)
           : (D == DOMAIN_ISO_TET)  ?((O == ORDER_CONSTANT       ) ? 1
                                    : (O == ORDER_LINEAR         ) ? 1
                                    : (O == ORDER_QUADRATIC      ) ? 4
                                    : (O == ORDER_CUBIC          ) ? 5
                                                                   : 0)
           : (D == DOMAIN_ISO_WEDGE)?((O == ORDER_SUPER_LINEAR   ) ? 2
                                    : (O == ORDER_SUPER_QUADRATIC) ? 9
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
            if (o == ORDER_LINEAR || o == ORDER_CONSTANT){
                points = {Point(1.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0)};
            }

            else if (o == ORDER_QUADRATIC){
                points[0] = Point(1.0/6.0, 1.0/6.0, 1.0/6.0),
                points[1] = Point(2.0/3.0, 1.0/6.0, 1.0/6.0),
                points[2] = Point(1.0/6.0, 2.0/3.0, 1.0/6.0);
            }

            else if (o == ORDER_CUBIC){
                points[0] = Point(1.0/3.0, 1.0/3.0, -27.0/96.0);
                points[1] = Point(1.0/5.0, 3.0/5.0,  25.0/96.0);
                points[2] = Point(1.0/5.0, 1.0/5.0,  25.0/96.0);
                points[3] = Point(3.0/5.0, 1.0/5.0,  25.0/96.0);
            }

            else {
                runtime_check(false, "not supported");
            }
        }

        if (d == DOMAIN_ISO_QUAD){

            if (o == ORDER_LINEAR || o == ORDER_CONSTANT){
                points[0] = Point(0.0, 0.0, 0.0, 4.0);
            }

            else if (o == ORDER_QUADRATIC || o == ORDER_CUBIC){
                points[0] = Point( sqrt(3.0)/3.0,  sqrt(3.0)/3.0, 1.0);
                points[1] = Point(-sqrt(3.0)/3.0,  sqrt(3.0)/3.0, 1.0);
                points[2] = Point(-sqrt(3.0)/3.0, -sqrt(3.0)/3.0, 1.0);
                points[3] = Point( sqrt(3.0)/3.0, -sqrt(3.0)/3.0, 1.0);
            }

            else if (o == ORDER_QUARTIC || o == ORDER_QUINTIC){
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

            else {
                runtime_check(false, "not supported");
            }
        }

        if (d == DOMAIN_ISO_TET) {
            if (o == ORDER_LINEAR || o == ORDER_CONSTANT){
                points[0] = Point(1 / 4.0, 1 / 4.0, 1 / 4.0, 1.0 / 6.0);
            }

            else if (o == ORDER_QUADRATIC) {
                Precision p1 = 1 / (3 * sqrt(5) - 5);
                Precision p2 = 1 / (    sqrt(5) + 5);
                points[0] = Point(p2, p2, p2, 1.0 / 4.0 / 6.0);
                points[1] = Point(p1, p2, p2, 1.0 / 4.0 / 6.0);
                points[2] = Point(p2, p1, p2, 1.0 / 4.0 / 6.0);
                points[3] = Point(p2, p2, p1, 1.0 / 4.0 / 6.0);
            }
            
            else if (o == ORDER_CUBIC) {
                points[0] = Point(1/4.0, 1/4.0, 1/4.0, -4/5.0 / 6);
                points[1] = Point(1/2.0, 1/6.0, 1/6.0, 9/20.0 / 6);
                points[2] = Point(1/6.0, 1/6.0, 1/6.0, 9/20.0 / 6);
                points[3] = Point(1/6.0, 1/6.0, 1/2.0, 9/20.0 / 6);
                points[4] = Point(1/6.0, 1/2.0, 1/6.0, 9/20.0 / 6);
            }

            else {
                runtime_check(false, "not supported");
            }
        }

        if (d == DOMAIN_ISO_HEX){
            if (o == ORDER_LINEAR || o == ORDER_CONSTANT){
                points[0] = Point(0.0, 0.0, 0.0, 8.0);
            }

            else if (o == ORDER_QUADRATIC || o == ORDER_CUBIC){
                points[0] = Point( 1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0), 1.0);
                points[1] = Point( 1.0/sqrt(3.0),  1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0);
                points[2] = Point( 1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0), 1.0);
                points[3] = Point( 1.0/sqrt(3.0), -1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0);
                points[4] = Point(-1.0/sqrt(3.0),  1.0/sqrt(3.0),  1.0/sqrt(3.0), 1.0);
                points[5] = Point(-1.0/sqrt(3.0),  1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0);
                points[6] = Point(-1.0/sqrt(3.0), -1.0/sqrt(3.0),  1.0/sqrt(3.0), 1.0);
                points[7] = Point(-1.0/sqrt(3.0), -1.0/sqrt(3.0), -1.0/sqrt(3.0), 1.0);
            }

            else if (o == ORDER_QUARTIC || o == ORDER_QUINTIC){
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

            else {
                runtime_check(false, "not supported");
            }
        }

        if (d == DOMAIN_ISO_WEDGE) {
            if(o == ORDER_SUPER_LINEAR){
                points[0] = Point(1.0 / 3.0, 1.0 / 3.0, - sqrt(3.0)/3.0 ,1.0 / 2.0);
                points[1] = Point(1.0 / 3.0, 1.0 / 3.0,   sqrt(3.0)/3.0 ,1.0 / 2.0);
            }

            else if(o == ORDER_SUPER_QUADRATIC){
                points[0] = Point(1.0/6.0, 1.0/6.0, -sqrt(0.6), 1.0 / 6.0 * 5.0 / 9.0),
                points[1] = Point(2.0/3.0, 1.0/6.0, -sqrt(0.6), 1.0 / 6.0 * 5.0 / 9.0),
                points[2] = Point(1.0/6.0, 2.0/3.0, -sqrt(0.6), 1.0 / 6.0 * 5.0 / 9.0);

                points[3] = Point(1.0/6.0, 1.0/6.0,  sqrt(0.0), 1.0 / 6.0 * 8.0 / 9.0),
                points[4] = Point(2.0/3.0, 1.0/6.0,  sqrt(0.0), 1.0 / 6.0 * 8.0 / 9.0),
                points[5] = Point(1.0/6.0, 2.0/3.0,  sqrt(0.0), 1.0 / 6.0 * 8.0 / 9.0);

                points[6] = Point(1.0/6.0, 1.0/6.0,  sqrt(0.6), 1.0 / 6.0 * 5.0 / 9.0),
                points[7] = Point(2.0/3.0, 1.0/6.0,  sqrt(0.6), 1.0 / 6.0 * 5.0 / 9.0),
                points[8] = Point(1.0/6.0, 2.0/3.0,  sqrt(0.6), 1.0 / 6.0 * 5.0 / 9.0);

            }

            else {
                runtime_check(false, "not supported");
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
