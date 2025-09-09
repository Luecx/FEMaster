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
        Precision It = get_profile()->I_t;
        Precision L = length();

        Precision a = E * A / L;
        Precision b = G * It / L;
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

    StaticMatrix<12,12> stiffness_geom_impl(IPData& ip_stress, int offset) {
        throw std::runtime_error("Not yet implemented");
    }


    StaticMatrix<12,12> mass_impl() override {
        StaticMatrix<12,12> T = transformation();
        StaticMatrix<12,12> M = StaticMatrix<12,12>::Zero();

        Precision A = get_profile()->A;
        Precision L = length();
        Precision rho = get_material()->get_density();

        Precision Ip = get_profile()->I_y + get_profile()->I_z;
        Precision IpA = Ip / A;

        // New stiffness matrix
        M <<
            140,        0,         0,         0,         0,         0,        70,         0,         0,         0,          0,          0,
              0,      156,         0,         0,         0,    22 * L,         0,        54,         0,         0,          0,    -13 * L,
              0,        0,       156,         0,   -22 * L,         0,         0,         0,        54,         0,     13 * L,          0,
              0,        0,         0, 140 * IpA,         0,         0,         0,         0,         0,  70 * IpA,          0,          0,
              0,        0,   -22 * L,         0, 4 * L * L,         0,         0,         0,   -13 * L,         0, -3 * L * L,          0,
              0,   22 * L,         0,         0,         0, 4 * L * L,         0,    13 * L,         0,         0,          0, -3 * L * L,
             70,        0,         0,         0,         0,         0,       140,         0,         0,         0,          0,          0,
              0,       54,         0,         0,         0,    13 * L,         0,       156,         0,         0,          0,    -22 * L,
              0,        0,        54,         0,   -13 * L,         0,         0,         0,       156,         0,     22 * L,          0,
              0,        0,         0,  70 * IpA,         0,         0,         0,         0,         0, 140 * IpA,          0,          0,
              0,        0,    13 * L,         0,-3 * L * L,         0,         0,         0,    22 * L,         0,  4 * L * L,          0,
              0,  -13 * L,         0,         0,         0,-3 * L * L,         0,   -22 * L,         0,         0,          0,  4 * L * L;

        M *= rho * L  * A / 420;

        return T.transpose() * M * T;
    }




};

}