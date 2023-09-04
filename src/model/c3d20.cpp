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
    Precision rp = r + 1;
    Precision sp = s + 1;
    Precision tp = t + 1;
    Precision rm = r - 1;
    Precision sm = s - 1;
    Precision tm = t - 1;

    // containing derivatives of the shape functions h1 -> h8 with respect to r/s/t
    StaticMatrix<20, 3> local_shape_derivative {};
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

        local_shape_derivative(i, 0) =
            -(Precision)1.0 / (Precision)8.0 * sum_s * sum_t * (sign_r1 * sum_1 + sign_r2 * sum_r);
        local_shape_derivative(i, 1) =
            -(Precision)1.0 / (Precision)8.0 * sum_r * sum_t * (sign_s1 * sum_1 + sign_s2 * sum_s);
        local_shape_derivative(i, 2) =
            -(Precision)1.0 / (Precision)8.0 * sum_r * sum_s * (sign_t1 * sum_1 + sign_t2 * sum_t);
    }

    // remaining 12 nodes
    for (int i = 8; i < 20; i++) {
        Precision r1      = (((i + 1) / 2) % 2 == 0 ? rm : rp);
        Precision r2      = (((i % 2 == 0 && i < 16 ? (i % 4 == 0 ? rp : rm) : 1)));
        Precision s1      = (((i) / 2) % 2 == 0 ? sm : sp);
        Precision s2      = (((i % 2 == 1 && i < 16 ? (i % 4 == 1 ? sp : sm) : 1)));
        Precision t1      = (((i) / 4) % 2 == 0 ? tm : tp);
        Precision t2      = (i < 16 ? 1 : tp);

        Precision sign_r1 = (Precision)(((i + 1) / 2) % 2 == 0 ? -1 : +1);
        Precision sign_r2 = (Precision)(((i % 2 == 0 && i < 16 ? (i % 4 == 0 ? +1 : -1) : 0)));
        Precision sign_s1 = (Precision)(((i) / 2) % 2 == 0 ? -1 : +1);
        Precision sign_s2 = (Precision)(((i % 2 == 1 && i < 16 ? (i % 4 == 1 ? +1 : -1) : 0)));
        Precision sign_t1 = (Precision)(((i) / 4) % 2 == 0 ? -1 : +1);
        Precision sign_t2 = (Precision)(i < 16 ? 0 : +1);

        local_shape_derivative(i, 0) = (Precision)1.0 / (Precision)4.0 * s1 * s2 * t1 * t2 * (r1 * sign_r2 + r2 * sign_r1);
        local_shape_derivative(i, 1) = (Precision)1.0 / (Precision)4.0 * r1 * r2 * t1 * t2 * (s1 * sign_s2 + s2 * sign_s1);
        local_shape_derivative(i, 2) = (Precision)1.0 / (Precision)4.0 * r1 * r2 * s1 * s2 * (t1 * sign_t2 + t2 * sign_t1);
    }

    return local_shape_derivative;
}

StaticMatrix<20, 3> C3D20::node_coords_local() {
    StaticMatrix<20, 3> res {};
    for (int n = 0; n < 8; n++) {
        Precision r1 = ((n + 1) / 2) % 2 == 0 ? -1 : 1;
        Precision s1 = ((n) / 2) % 2 == 0 ? -1 : 1;
        Precision t1 = n >= 4 ? 1 : -1;

        res(n, 0)    = r1;
        res(n, 1)    = s1;
        res(n, 2)    = t1;
    }
    return res;
}

StaticMatrix<6, 60> C3D20::strain_displacement(const StaticMatrix<20, 3>& shape_der_global) {
    StaticMatrix<6, 60> B {};
    B.setZero();
    for (int j = 0; j < 20; j++) {
        int r1   = j * 3;
        int r2   = r1 + 1;
        int r3   = r1 + 2;
        B(0, r1) = shape_der_global(j, 0);
        B(1, r2) = shape_der_global(j, 1);
        B(2, r3) = shape_der_global(j, 2);
        B(3, r1) = shape_der_global(j, 1) / (Precision) 2.0;    // divide by 2 in order to account for real shear strain
        B(3, r2) = shape_der_global(j, 0) / (Precision) 2.0;    // divide by 2 in order to account for real shear strain
        B(4, r1) = shape_der_global(j, 2) / (Precision) 2.0;    // divide by 2 in order to account for real shear strain
        B(4, r3) = shape_der_global(j, 0) / (Precision) 2.0;    // divide by 2 in order to account for real shear strain
        B(5, r2) = shape_der_global(j, 2) / (Precision) 2.0;    // divide by 2 in order to account for real shear strain
        B(5, r3) = shape_der_global(j, 1) / (Precision) 2.0;    // divide by 2 in order to account for real shear strain
    }
    return B;
}

}    // namespace model
}    // namespace fem