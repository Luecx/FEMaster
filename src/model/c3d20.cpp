//
// Created by Luecx on 12.06.2023.
//

#include "c3d20.h"

namespace fem {
namespace model {

C3D20::C3D20(ID p_elem_id, const std::array<ID, 20>& p_node_ids)
    : SolidElement(p_elem_id, p_node_ids) {}

StaticMatrix<20, 1> C3D20::shape_function(Precision r, Precision s, Precision t) {
    StaticMatrix<20, 1>  res {};

    Precision          rp = 1 + r;
    Precision          sp = 1 + s;
    Precision          tp = 1 + t;
    Precision          rm = 1 - r;
    Precision          sm = 1 - s;
    Precision          tm = 1 - t;

    // first 8 nodes
    for (int i = 0; i < 8; i++) {
        Precision sign_r1 = (Precision) (((i + 1) / 2) % 2 == 0 ? -1 : 1);
        Precision sign_s1 = (Precision) ((i / 2) % 2 == 0 ? -1 : 1);
        Precision sign_t1 = (Precision) ((i / 4) % 2 == 0 ? -1 : 1);

        Precision sign_r2 = -sign_r1;
        Precision sign_s2 = -sign_s1;
        Precision sign_t2 = -sign_t1;

        Precision sum_1   = 2 + r * sign_r2 + s * sign_s2 + t * sign_t2;

        Precision sum_r   = (((i + 1) / 2) % 2 == 0 ? rm : rp);
        Precision sum_s   = ((i / 2) % 2 == 0 ? sm : sp);
        Precision sum_t   = ((i / 4) % 2 == 0 ? tm : tp);

        res(i, 0)         = (Precision) -1.0 / (Precision) 8.0 * sum_r * sum_s * sum_t * sum_1;
    }
    // remaining 12 nodes
    for (int i = 8; i < 20; i++) {
        Precision r1 = (((i + 1) / 2) % 2 == 0 ? rm : rp);
        Precision r2 = (((i % 2 == 0 && i < 16 ? (i % 4 == 0 ? rp : rm) : 1)));
        Precision s1 = (((i) / 2) % 2 == 0 ? sm : sp);
        Precision s2 = (((i % 2 == 1 && i < 16 ? (i % 4 == 1 ? sp : sm) : 1)));
        Precision t1 = (((i) / 4) % 2 == 0 ? tm : tp);
        Precision t2 = (i < 16 ? 1 : tp);

        res(i, 0)    = (Precision) 1.0 / (Precision) 4.0 * r1 * r2 * s1 * s2 * t1 * t2;
    }
    return res;
}

StaticMatrix<20, 3> C3D20::shape_derivative(Precision r, Precision s, Precision t) {

    Precision g = r;
    Precision h = s;
    r = t;

    Precision gp = 1 + g;
    Precision hp = 1 + h;
    Precision rp = 1 + r;
    Precision gm = 1 - g;
    Precision hm = 1 - h;
    Precision rm = 1 - r;

    // containing derivatives of the shape functions h1 -> h8 with respect to r/s/t
    StaticMatrix<20, 3> der {};
    der(0,0) = - 1.0 / 8.0 * (hm * rm) * (-1 * (2 + g + h + r) + 1 * gm);
    der(1,0) = - 1.0 / 8.0 * (hm * rm) * ( 1 * (2 - g + h + r) - 1 * gp);
    der(2,0) = - 1.0 / 8.0 * (hp * rm) * ( 1 * (2 - g - h + r) - 1 * gp);
    der(3,0) = - 1.0 / 8.0 * (hp * rm) * (-1 * (2 + g - h + r) + 1 * gm);
    der(4,0) = - 1.0 / 8.0 * (hm * rp) * (-1 * (2 + g + h - r) + 1 * gm);
    der(5,0) = - 1.0 / 8.0 * (hm * rp) * ( 1 * (2 - g + h - r) - 1 * gp);
    der(6,0) = - 1.0 / 8.0 * (hp * rp) * ( 1 * (2 - g - h - r) - 1 * gp);
    der(7,0) = - 1.0 / 8.0 * (hp * rp) * (-1 * (2 + g - h - r) + 1 * gm);

    der(0,1) = - 1.0 / 8.0 * (gm * rm) * (-1 * (2 + g + h + r) + 1 * hm);
    der(1,1) = - 1.0 / 8.0 * (gp * rm) * (-1 * (2 - g + h + r) + 1 * hm);
    der(2,1) = - 1.0 / 8.0 * (gp * rm) * ( 1 * (2 - g - h + r) - 1 * hp);
    der(3,1) = - 1.0 / 8.0 * (gm * rm) * ( 1 * (2 + g - h + r) - 1 * hp);
    der(4,1) = - 1.0 / 8.0 * (gm * rp) * (-1 * (2 + g + h - r) + 1 * hm);
    der(5,1) = - 1.0 / 8.0 * (gp * rp) * (-1 * (2 - g + h - r) + 1 * hm);
    der(6,1) = - 1.0 / 8.0 * (gp * rp) * ( 1 * (2 - g - h - r) - 1 * hp);
    der(7,1) = - 1.0 / 8.0 * (gm * rp) * ( 1 * (2 + g - h - r) - 1 * hp);

    der(0,2) = - 1.0 / 8.0 * (gm * hm) * (-1 * (2 + g + h + r) + 1 * rm);
    der(1,2) = - 1.0 / 8.0 * (gp * hm) * (-1 * (2 - g + h + r) + 1 * rm);
    der(2,2) = - 1.0 / 8.0 * (gp * hp) * (-1 * (2 - g - h + r) + 1 * rm);
    der(3,2) = - 1.0 / 8.0 * (gm * hp) * (-1 * (2 + g - h + r) + 1 * rm);
    der(4,2) = - 1.0 / 8.0 * (gm * hm) * ( 1 * (2 + g + h - r) - 1 * rp);
    der(5,2) = - 1.0 / 8.0 * (gp * hm) * ( 1 * (2 - g + h - r) - 1 * rp);
    der(6,2) = - 1.0 / 8.0 * (gp * hp) * ( 1 * (2 - g - h - r) - 1 * rp);
    der(7,2) = - 1.0 / 8.0 * (gm * hp) * ( 1 * (2 + g - h - r) - 1 * rp);

    // those were g comes twice
    der(8 ,0) = 1.0 / 4.0 * (hm * rm) * ( 1 * gm - 1 * gp );
    der(10,0) = 1.0 / 4.0 * (hp * rm) * ( 1 * gm - 1 * gp );
    der(12,0) = 1.0 / 4.0 * (hm * rp) * ( 1 * gm - 1 * gp );
    der(14,0) = 1.0 / 4.0 * (hp * rp) * ( 1 * gm - 1 * gp );

    // those were h comes twice
    der(9 ,1) = 1.0 / 4.0 * (gp * rm) * ( 1 * hm - 1 * hp );
    der(11,1) = 1.0 / 4.0 * (gm * rm) * ( 1 * hm - 1 * hp );
    der(13,1) = 1.0 / 4.0 * (gp * rp) * ( 1 * hm - 1 * hp );
    der(15,1) = 1.0 / 4.0 * (gm * rp) * ( 1 * hm - 1 * hp );

    // those were r comes twice
    der(16,2) = 1.0 / 4.0 * (gm * hm) * ( 1 * rm - 1 * rp );
    der(17,2) = 1.0 / 4.0 * (gp * hm) * ( 1 * rm - 1 * rp );
    der(18,2) = 1.0 / 4.0 * (gp * hp) * ( 1 * rm - 1 * rp );
    der(19,2) = 1.0 / 4.0 * (gm * hp) * ( 1 * rm - 1 * rp );


    // those were g comes twice
    der(8 ,1) = 1.0 / 4.0 * (-1 * rm) * (gm * gp);
    der(10,1) = 1.0 / 4.0 * ( 1 * rm) * (gm * gp);
    der(12,1) = 1.0 / 4.0 * (-1 * rp) * (gm * gp);
    der(14,1) = 1.0 / 4.0 * ( 1 * rp) * (gm * gp);
    der(8 ,2) = 1.0 / 4.0 * (-1 * hm) * (gm * gp);
    der(10,2) = 1.0 / 4.0 * (-1 * hp) * (gm * gp);
    der(12,2) = 1.0 / 4.0 * ( 1 * hm) * (gm * gp);
    der(14,2) = 1.0 / 4.0 * ( 1 * hp) * (gm * gp);

    // those were h comes twice
    der(9 ,0) = 1.0 / 4.0 * ( 1 * rm) * (hm * hp);
    der(11,0) = 1.0 / 4.0 * (-1 * rm) * (hm * hp);
    der(13,0) = 1.0 / 4.0 * ( 1 * rp) * (hm * hp);
    der(15,0) = 1.0 / 4.0 * (-1 * rp) * (hm * hp);
    der(9 ,2) = 1.0 / 4.0 * (-1 * gp) * (hm * hp);
    der(11,2) = 1.0 / 4.0 * (-1 * gm) * (hm * hp);
    der(13,2) = 1.0 / 4.0 * ( 1 * gp) * (hm * hp);
    der(15,2) = 1.0 / 4.0 * ( 1 * gm) * (hm * hp);

    // those were r comes twice
    der(16,0) = 1.0 / 4.0 * (-1 * hm) * (rm * rp);
    der(17,0) = 1.0 / 4.0 * ( 1 * hm) * (rm * rp);
    der(18,0) = 1.0 / 4.0 * ( 1 * hp) * (rm * rp);
    der(19,0) = 1.0 / 4.0 * (-1 * hp) * (rm * rp);
    der(16,1) = 1.0 / 4.0 * (-1 * gm) * (rm * rp);
    der(17,1) = 1.0 / 4.0 * (-1 * gp) * (rm * rp);
    der(18,1) = 1.0 / 4.0 * ( 1 * gp) * (rm * rp);
    der(19,1) = 1.0 / 4.0 * ( 1 * gm) * (rm * rp);


//    // first 8 nodes
//    for (int i = 0; i < 8; i++) {
//        Precision sign_r1 = (Precision) (((i + 1) / 2) % 2 == 0 ? -1 : 1);
//        Precision sign_s1 = (Precision) ((i / 2) % 2 == 0 ? -1 : 1);
//        Precision sign_t1 = (Precision) ((i / 4) % 2 == 0 ? -1 : 1);
//
//        Precision sign_r2 = -sign_r1;
//        Precision sign_s2 = -sign_s1;
//        Precision sign_t2 = -sign_t1;
//
//        Precision sum_1   = 2 + r * sign_r2 + s * sign_s2 + t * sign_t2;
//
//        Precision sum_r   = (((i + 1) / 2) % 2 == 0 ? rm : rp);
//        Precision sum_s   = ((i / 2) % 2 == 0 ? sm : sp);
//        Precision sum_t   = ((i / 4) % 2 == 0 ? tm : tp);
//
//        local_shape_derivative(i, 0) =
//            -(Precision)1.0 / (Precision)8.0 * sum_s * sum_t * (sign_r1 * sum_1 + sign_r2 * sum_r);
//        local_shape_derivative(i, 1) =
//            -(Precision)1.0 / (Precision)8.0 * sum_r * sum_t * (sign_s1 * sum_1 + sign_s2 * sum_s);
//        local_shape_derivative(i, 2) =
//            -(Precision)1.0 / (Precision)8.0 * sum_r * sum_s * (sign_t1 * sum_1 + sign_t2 * sum_t);
//    }
//
//    // remaining 12 nodes
//    for (int i = 8; i < 20; i++) {
//        Precision r1      = (((i + 1) / 2) % 2 == 0 ? rm : rp);
//        Precision r2      = (((i % 2 == 0 && i < 16 ? (i % 4 == 0 ? rp : rm) : 1)));
//        Precision s1      = (((i) / 2) % 2 == 0 ? sm : sp);
//        Precision s2      = (((i % 2 == 1 && i < 16 ? (i % 4 == 1 ? sp : sm) : 1)));
//        Precision t1      = (((i) / 4) % 2 == 0 ? tm : tp);
//        Precision t2      = (i < 16 ? 1 : tp);
//
//        Precision sign_r1 = (Precision)(((i + 1) / 2) % 2 == 0 ? -1 : +1);
//        Precision sign_r2 = (Precision)(((i % 2 == 0 && i < 16 ? (i % 4 == 0 ? +1 : -1) : 0)));
//        Precision sign_s1 = (Precision)(((i) / 2) % 2 == 0 ? -1 : +1);
//        Precision sign_s2 = (Precision)(((i % 2 == 1 && i < 16 ? (i % 4 == 1 ? +1 : -1) : 0)));
//        Precision sign_t1 = (Precision)(((i) / 4) % 2 == 0 ? -1 : +1);
//        Precision sign_t2 = (Precision)(i < 16 ? 0 : +1);
//
//        local_shape_derivative(i, 0) = (Precision)1.0 / (Precision)4.0 * s1 * s2 * t1 * t2 * (r1 * sign_r2 + r2 * sign_r1);
//        local_shape_derivative(i, 1) = (Precision)1.0 / (Precision)4.0 * r1 * r2 * t1 * t2 * (s1 * sign_s2 + s2 * sign_s1);
//        local_shape_derivative(i, 2) = (Precision)1.0 / (Precision)4.0 * r1 * r2 * s1 * s2 * (t1 * sign_t2 + t2 * sign_t1);
//    }

    return der;
}

StaticMatrix<20, 3> C3D20::node_coords_local() {
    StaticMatrix<20, 3> res {};
    res.setZero();
    for (int n = 0; n < 8; n++) {
        Precision r1 = ((n + 1) / 2) % 2 == 0 ? -1 : 1;
        Precision s1 = ((n) / 2) % 2 == 0 ? -1 : 1;
        Precision t1 = n >= 4 ? 1 : -1;

        res(n, 0)    = r1;
        res(n, 1)    = s1;
        res(n, 2)    = t1;
    }

    // Mid-edge nodes
    res( 8, 0) =  0.0; res( 8, 1) = -1.0; res( 8, 2) = -1.0;
    res( 9, 0) =  1.0; res( 9, 1) =  0.0; res( 9, 2) = -1.0;
    res(10, 0) =  0.0; res(10, 1) =  1.0; res(10, 2) = -1.0;
    res(11, 0) = -1.0; res(11, 1) =  0.0; res(11, 2) = -1.0;
    res(12, 0) =  0.0; res(12, 1) = -1.0; res(12, 2) =  1.0;
    res(13, 0) =  1.0; res(13, 1) =  0.0; res(13, 2) =  1.0;
    res(14, 0) =  0.0; res(14, 1) =  1.0; res(14, 2) =  1.0;
    res(15, 0) = -1.0; res(15, 1) =  0.0; res(15, 2) =  1.0;
    res(16, 0) = -1.0; res(16, 1) = -1.0; res(16, 2) =  0.0;
    res(17, 0) =  1.0; res(17, 1) = -1.0; res(17, 2) =  0.0;
    res(18, 0) =  1.0; res(18, 1) =  1.0; res(18, 2) =  0.0;
    res(19, 0) = -1.0; res(19, 1) =  1.0; res(19, 2) =  0.0;

    return res;
}

}    // namespace model
}    // namespace fem