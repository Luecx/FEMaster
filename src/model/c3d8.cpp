//
// Created by Luecx on 12.06.2023.
//

#include "c3d8.h"

namespace fem {
namespace model {

C3D8::C3D8(ID pElemId, const std::array<ID, 8>& pNodeIds)
    : SolidElement(pElemId, pNodeIds) {}

StaticMatrix<8, 1> C3D8::shape_function(Precision r, Precision s, Precision t) {
    StaticMatrix<8, 1>  res {};

    Precision rp = r + 1;
    Precision sp = s + 1;
    Precision tp = t + 1;
    Precision rm = r - 1;
    Precision sm = s - 1;
    Precision tm = t - 1;

    // shape function evaluation. 4 entries for 1) sign, 2) r 3) s 4) t
    for (int n = 0; n < 8; n++) {
        Precision sign = (n == 0 || n == 2 || n == 5 || n == 7) ? -1 : 1;
        Precision r1   = ((n + 1) / 2) % 2 == 1 ? rp : rm;
        Precision s1   = ((n) / 2) % 2 == 1 ? sp : sm;
        Precision t1   = n >= 4 ? tp : tm;

        res(n, 0)      = sign * r1 * s1 * t1;
    }

    res *= 0.125;
    return res;
}

StaticMatrix<8, 3> C3D8::shape_derivative(Precision r, Precision s, Precision t) {
    Precision rp = r + 1;
    Precision sp = s + 1;
    Precision tp = t + 1;
    Precision rm = r - 1;
    Precision sm = s - 1;
    Precision tm = t - 1;

    // shape function evaluation. 4 entries for 1) sign, 2) r 3) s 4) t
    StaticMatrix<8, 4> shape_function {};
    for (int n = 0; n < 8; n++) {
        shape_function(n, 0) = (n == 0 || n == 2 || n == 5 || n == 7) ? -1 : 1;
        shape_function(n, 1) = ((n + 1) / 2) % 2 == 1 ? rp : rm;
        shape_function(n, 2) = ((n) / 2) % 2 == 1 ? sp : sm;
        shape_function(n, 3) = n >= 4 ? tp : tm;
    }

    // containing derivatives of the shape functions h1 -> h8 with respect to r/s/t
    StaticMatrix<8, 3> local_shape_derivative {};
    for (int n = 0; n < 8; n++) {
        local_shape_derivative(n, 0) = shape_function(n, 0) * shape_function(n, 2) * shape_function(n, 3);
        local_shape_derivative(n, 1) = shape_function(n, 0) * shape_function(n, 1) * shape_function(n, 3);
        local_shape_derivative(n, 2) = shape_function(n, 0) * shape_function(n, 1) * shape_function(n, 2);
    }
    local_shape_derivative *= 0.125;
    return local_shape_derivative;
}

StaticMatrix<8, 3> C3D8::node_coords_local() {
    StaticMatrix<8, 3> res {};
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

StaticMatrix<6, 24> C3D8::strain_displacement(const StaticMatrix<8, 3>& shape_der_global) {
    StaticMatrix<6, 24> B {};
    B.setZero();
    for (int j = 0; j < 8; j++) {
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