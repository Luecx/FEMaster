#pragma once

#include "beam.h"

namespace fem::model {

struct B33 : BeamElement<2>{

    B33 (ID p_elem_id, std::array<ID, 2> p_node_ids)
        : BeamElement(p_elem_id, p_node_ids) {}

    StaticMatrix<12,12> stiffness_impl() override {

        StaticMatrix<12,12> T = transformation();
        StaticMatrix<12,12> K = StaticMatrix<12,12>::Zero();

        Precision E = get_elasticity()->youngs;
        Precision G = get_elasticity()->shear;
        Precision A = get_profile()->A;
        Precision Iy = get_profile()->I_y;
        Precision Iz = get_profile()->I_z;
        Precision Jt = get_profile()->J_t;
        Precision L = length();

        Precision a = E * A / L;
        Precision b = G * Jt / L;
        Precision c = E * Iz / (L * L * L);
        Precision d = E * Iy / (L * L * L);

        Precision M = L * L;

        K <<
            a,         0,         0,        0,         0,         0,       -a,         0,         0,        0,         0,         0,
            0,    12 * c,         0,        0,         0, 6 * c * L,        0,   -12 * c,         0,        0,         0, 6 * c * L,
            0,         0,    12 * d,        0,-6 * d * L,         0,        0,         0,   -12 * d,        0,-6 * d * L,         0,
            0,         0,         0,        b,         0,         0,        0,         0,         0,       -b,         0,         0,
            0,         0,-6 * d * L,        0, 4 * d * M,         0,        0,         0, 6 * d * L,        0, 2 * d * M,         0,
            0, 6 * c * L,         0,        0,         0, 4 * c * M,        0,-6 * c * L,         0,        0,         0, 2 * c * M,
           -a,         0,         0,        0,         0,         0,        a,         0,         0,        0,         0,         0,
            0,   -12 * c,         0,        0,         0,-6 * c * L,        0,    12 * c,         0,        0,         0,-6 * c * L,
            0,         0,   -12 * d,        0, 6 * d * L,         0,        0,         0,    12 * d,        0, 6 * d * L,         0,
            0,         0,         0,       -b,         0,         0,        0,         0,         0,        b,         0,         0,
            0,         0,-6 * d * L,        0, 2 * d * M,         0,        0,         0, 6 * d * L,        0, 4 * d * M,         0,
            0, 6 * c * L,         0,        0,         0, 2 * c * M,        0,-6 * c * L,         0,        0,         0, 4 * c *M;

        return T.transpose() * K * T;
    }

    StaticMatrix<12,12> mass_impl() override {
        return StaticMatrix<12,12>::Zero();
    }




};

}