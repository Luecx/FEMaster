//
// Created by Luecx on 30.09.2023.
//

#include "quadrature.h"
#include "../core/assert.h"
#include <iostream>


namespace fem{

namespace quadrature {

constexpr size_t n_dims(Order O, Domain D) {
    return D == DOMAIN_ISO_TRI ? 2
         : D == DOMAIN_ISO_QUAD ? 2
         : D == DOMAIN_ISO_HEX ? 3
         : D == DOMAIN_ISO_TET ? 3
         : D == DOMAIN_ISO_WEDGE ? 3
         : D == DOMAIN_ISO_LINE_A ? 1
         : D == DOMAIN_ISO_LINE_B ? 1
         : D == DOMAIN_ISO_PYRAMID ? 3
         : 3;
}
constexpr size_t n_points(Order O, Domain D) {
    return   (D == DOMAIN_ISO_TRI)       ? ((O == ORDER_CONSTANT       ) ? 1
                                          : (O == ORDER_LINEAR         ) ? 1
                                          : (O == ORDER_QUADRATIC      ) ? 3
                                          : (O == ORDER_CUBIC          ) ? 4
                                                                         : 0)
           : (D == DOMAIN_ISO_QUAD)      ? ((O == ORDER_CONSTANT       ) ? 1
                                          : (O == ORDER_LINEAR         ) ? 1
                                          : (O == ORDER_QUADRATIC      ) ? 4
                                          : (O == ORDER_CUBIC          ) ? 4
                                          : (O == ORDER_QUARTIC        ) ? 9
                                          : (O == ORDER_QUINTIC        ) ? 9
                                                                         : 0)
           : (D == DOMAIN_ISO_HEX)       ? ((O == ORDER_CONSTANT       ) ? 1
                                          : (O == ORDER_LINEAR         ) ? 1
                                          : (O == ORDER_QUADRATIC      ) ? 8
                                          : (O == ORDER_CUBIC          ) ? 8
                                          : (O == ORDER_QUARTIC        ) ? 27
                                          : (O == ORDER_QUINTIC        ) ? 27
                                                                         : 0)
           : (D == DOMAIN_ISO_PYRAMID)   ? ((O == ORDER_CONSTANT       ) ? 1
                                          : (O == ORDER_LINEAR         ) ? 1
                                          : (O == ORDER_QUADRATIC      ) ? 8
                                          : (O == ORDER_CUBIC          ) ? 8
                                          : (O == ORDER_QUARTIC        ) ? 27
                                          : (O == ORDER_QUINTIC        ) ? 27
                                                                         : 0)
           : (D == DOMAIN_ISO_TET)     ? ((O == ORDER_CONSTANT       ) ? 1
                                          : (O == ORDER_LINEAR         ) ? 1
                                          : (O == ORDER_QUADRATIC      ) ? 4
                                          : (O == ORDER_CUBIC          ) ? 5
                                          : (O == ORDER_QUARTIC        ) ? 8
                                          : (O == ORDER_QUINTIC        ) ? 15
                                                                         : 0)
           : (D == DOMAIN_ISO_WEDGE)     ? ((O == ORDER_SUPER_LINEAR   ) ? 2
                                          : (O == ORDER_SUPER_QUADRATIC) ? 9
                                                                         : 0)
           : (D == DOMAIN_ISO_LINE_A)    ? ((O == ORDER_CONSTANT       ) ? 1
                                          : (O == ORDER_LINEAR         ) ? 1
                                          : (O == ORDER_QUADRATIC      ) ? 2
                                          : (O == ORDER_CUBIC          ) ? 2
                                          : (O == ORDER_QUARTIC        ) ? 3
                                          : (O == ORDER_QUINTIC        ) ? 3
                                                                         : 0)
           : (D == DOMAIN_ISO_LINE_B)    ? ((O == ORDER_CONSTANT       ) ? 1
                                          : (O == ORDER_LINEAR         ) ? 1
                                          : (O == ORDER_QUADRATIC      ) ? 2
                                          : (O == ORDER_CUBIC          ) ? 2
                                          : (O == ORDER_QUARTIC        ) ? 3
                                          : (O == ORDER_QUINTIC        ) ? 3
                                                                         : 0)
                                       : 0;
}

Point::Point(Precision r, Precision s, Precision t, Precision w)
    : r(r)
    , s(s)
    , t(t)
    , w(w) {}
Point::Point(Precision r, Precision s, Precision w)
    : r(r)
    , s(s)
    , t(0)
    , w(w) {}
Point::Point(Precision r, Precision w)
    : r(r)
    , s(0)
    , t(0)
    , w(w) {}

Quadrature::Quadrature(Domain d, Order o)
    : domain(d), order(o) {
    points.resize(n_points(o, d));

    if (d == DOMAIN_ISO_LINE_A) {
        if (o == ORDER_CONSTANT || o == ORDER_LINEAR)
            points = {Point(0.0, 2.0)};
        else if (o == ORDER_QUADRATIC || o == ORDER_CUBIC){
            points[0] = Point(-1.0/sqrt(3.0), 1.0);
            points[1] = Point( 1.0/sqrt(3.0), 1.0);
        }
        else if (o == ORDER_QUARTIC || o == ORDER_QUINTIC){
            points[0] = Point(-sqrt(3.0/5.0), 5.0/9.0);
            points[1] = Point( 0.0, 8.0/9.0);
            points[2] = Point( sqrt(3.0/5.0), 5.0/9.0);
        }
        else {
            runtime_check(false, "not supported");
        }
    }

    if (d == DOMAIN_ISO_LINE_B) {
        if (o == ORDER_CONSTANT || o == ORDER_LINEAR)
            points = {Point(0.5, 1.0)};
        else if (o == ORDER_QUADRATIC || o == ORDER_CUBIC){
            points[0] = Point(0.5 - 0.5/sqrt(3.0), 0.5);
            points[1] = Point(0.5 + 0.5/sqrt(3.0), 0.5);
        }
        else if (o == ORDER_QUARTIC || o == ORDER_QUINTIC){
            points[0] = Point(0.5 - 0.5*sqrt(3.0/5.0), 5.0/18.0);
            points[1] = Point(0.5, 8.0/18.0);
            points[2] = Point(0.5 + 0.5*sqrt(3.0/5.0), 5.0/18.0);
        }
        else {
            runtime_check(false, "not supported");
        }
    }

    if (d == DOMAIN_ISO_TRI){
        if (o == ORDER_LINEAR || o == ORDER_CONSTANT){
            points = {Point(1.0 / 3.0, 1.0 / 3.0, 1.0 / 2.0)};
        }

        else if (o == ORDER_QUADRATIC){
            points[0] = Point(1.0/6.0, 1.0/6.0, 1.0/6.0);
            points[1] = Point(2.0/3.0, 1.0/6.0, 1.0/6.0);
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

        // New: 8-point Gauss rule for Quartic order
        else if (o == ORDER_QUARTIC) {
            points[0] = Point(0.0158359099, 0.3280546970, 0.3280546970, 0.138527967 / 6);
            points[1] = Point(0.3280546970, 0.0158359099, 0.3280546970, 0.138527967 / 6);
            points[2] = Point(0.3280546970, 0.3280546970, 0.0158359099, 0.138527967 / 6);
            points[3] = Point(0.3280546970, 0.3280546970, 0.3280546970, 0.138527967 / 6);
            points[4] = Point(0.6791431780, 0.1069522740, 0.1069522740, 0.111472033 / 6);
            points[5] = Point(0.1069522740, 0.6791431780, 0.1069522740, 0.111472033 / 6);
            points[6] = Point(0.1069522740, 0.1069522740, 0.6791431780, 0.111472033 / 6);
            points[7] = Point(0.1069522740, 0.1069522740, 0.1069522740, 0.111472033 / 6);
        }

        // 15-point rule for Quintic order
        else if (o == ORDER_QUINTIC) {
            points[0]  = Point(0.25, 0.25, 0.25, 0.030283678097089);
            points[1]  = Point(0.333333333333333, 0.333333333333333, 0.333333333333333, 0.006026785714286);
            points[2]  = Point(0.0, 0.333333333333333, 0.333333333333333, 0.006026785714286);
            points[3]  = Point(0.333333333333333, 0.0, 0.333333333333333, 0.006026785714286);
            points[4]  = Point(0.333333333333333, 0.333333333333333, 0.0, 0.006026785714286);
            points[5]  = Point(0.090909090909091, 0.090909090909091, 0.090909090909091, 0.011645249086029);
            points[6]  = Point(0.727272727272727, 0.090909090909091, 0.090909090909091, 0.011645249086029);
            points[7]  = Point(0.090909090909091, 0.727272727272727, 0.090909090909091, 0.011645249086029);
            points[8]  = Point(0.090909090909091, 0.090909090909091, 0.727272727272727, 0.011645249086029);
            points[9]  = Point(0.433449846426336, 0.066550153573664, 0.066550153573664, 0.010949141561386);
            points[10] = Point(0.066550153573664, 0.433449846426336, 0.066550153573664, 0.010949141561386);
            points[11] = Point(0.066550153573664, 0.066550153573664, 0.433449846426336, 0.010949141561386);
            points[12] = Point(0.066550153573664, 0.433449846426336, 0.433449846426336, 0.010949141561386);
            points[13] = Point(0.433449846426336, 0.066550153573664, 0.433449846426336, 0.010949141561386);
            points[14] = Point(0.433449846426336, 0.433449846426336, 0.066550153573664, 0.010949141561386);
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

    if (d == DOMAIN_ISO_PYRAMID){
        if (o == ORDER_LINEAR || o == ORDER_CONSTANT){
            points[0] = Point(0.0, 0.0, 1.0 / 3.0, 4.0 / 3.0);
        }

        else if (o == ORDER_QUADRATIC || o == ORDER_CUBIC) {
            Precision z_lower = 0.122514822655441;
            Precision z_upper = 0.544151844011225;
            Precision w_lower = 0.100785882079825;
            Precision w_upper = 0.232547451253500;

            Precision x_upper = 0.263184055569713;
            Precision x_lower = 0.506616303349787;

            points[0] = Point( x_upper,  x_upper, z_upper, w_upper);
            points[1] = Point(-x_upper,  x_upper, z_upper, w_upper);
            points[2] = Point(-x_upper, -x_upper, z_upper, w_upper);
            points[3] = Point( x_upper, -x_upper, z_upper, w_upper);

            points[4] = Point( x_lower,  x_lower, z_lower, w_lower);
            points[5] = Point(-x_lower,  x_lower, z_lower, w_lower);
            points[6] = Point(-x_lower, -x_lower, z_lower, w_lower);
            points[7] = Point( x_lower, -x_lower, z_lower, w_lower);
        }


        else if (o == ORDER_QUARTIC || o == ORDER_QUINTIC) {
            // Gr, Gs, and Gt values as arrays
            Precision gr[] = { -0.228504305653967, -0.505808707853924, -0.718055741319888,
                               -0.228504305653967, -0.505808707853924, -0.718055741319888,
                               -0.228504305653967, -0.505808707853924, -0.718055741319888,
                                0.000000000000000,  0.000000000000000,  0.000000000000000,
                                0.000000000000000,  0.000000000000000,  0.000000000000000,
                                0.000000000000000,  0.000000000000000,  0.000000000000000,
                                0.228504305653967,  0.505808707853924,  0.718055741319888,
                                0.228504305653967,  0.505808707853924,  0.718055741319888,
                                0.228504305653967,  0.505808707853924,  0.718055741319888 };

            Precision gs[] = { -0.228504305653967, -0.505808707853924, -0.718055741319888,
                                0.000000000000000,  0.000000000000000,  0.000000000000000,
                                0.228504305653967,  0.505808707853924,  0.718055741319888,
                               -0.228504305653967, -0.505808707853924, -0.718055741319888,
                                0.000000000000000,  0.000000000000000,  0.000000000000000,
                                0.228504305653967,  0.505808707853924,  0.718055741319888,
                               -0.228504305653967, -0.505808707853924, -0.718055741319888,
                                0.000000000000000,  0.000000000000000,  0.000000000000000,
                                0.228504305653967,  0.505808707853924,  0.718055741319888 };

            Precision gt[] = {  0.705002209888498,  0.347003766038352,  0.072994024073150,
                                0.705002209888498,  0.347003766038352,  0.072994024073150,
                                0.705002209888498,  0.347003766038352,  0.072994024073150,
                                0.705002209888498,  0.347003766038352,  0.072994024073150,
                                0.705002209888498,  0.347003766038352,  0.072994024073150,
                                0.705002209888498,  0.347003766038352,  0.072994024073150,
                                0.705002209888498,  0.347003766038352,  0.072994024073150,
                                0.705002209888498,  0.347003766038352,  0.072994024073150,
                                0.705002209888498,  0.347003766038352,  0.072994024073150 };

            // Weights (from the table)
            Precision weights[] = { 0.009244044138451, 0.045137737425885, 0.048498876871879,
                                    0.014790470621521, 0.072220379881415, 0.0775982029950066,
                                    0.009244044138451, 0.045137737425885, 0.048498876871879,
                                    0.014790470621521, 0.072220379881415, 0.0775982029950066,
                                    0.023664752994434, 0.115552607810264, 0.124157124792009,
                                    0.014790470621521, 0.072220379881415, 0.0775982029950066,
                                    0.009244044138451, 0.045137737425885, 0.048498876871879,
                                    0.014790470621521, 0.072220379881415, 0.0775982029950066,
                                    0.009244044138451, 0.045137737425885, 0.048498876871879};

            // Loop to fill in the points
            for (int i = 0; i < 27; ++i) {
                points[i] = Point(gr[i], gs[i], gt[i], weights[i]);
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
Point Quadrature::get_point(ID n) const {
    return points[n];
}
ID Quadrature::count() const {
    return points.size();
}
}
}