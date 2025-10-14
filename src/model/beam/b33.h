/**
 * @file b33.h
 * @brief Declares the two-node 3D beam element (B33).
 *
 * The element builds on `BeamElement<2>` and provides stiffness, geometric
 * stiffness, and mass matrices together with basic stress recovery.
 *
 * @see src/model/beam/beam.h
 */

#pragma once

#include "beam.h"

#include <limits>

namespace fem {
namespace model {

struct B33 : BeamElement<2> {
    B33(ID elem_id, std::array<ID, 2> node_ids_in)
        : BeamElement(elem_id, node_ids_in) {}

    StaticMatrix<12, 12> stiffness_impl() override {
        StaticMatrix<12, 12> T = transformation();
        StaticMatrix<12, 12> K = StaticMatrix<12, 12>::Zero();

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
            0, 6 * c * L,         0,        0,         0, 2 * c * M,        0,-6 * c * L,         0,        0,         0, 4 * c * M;

        return T.transpose() * K * T;
    }

    StaticMatrix<12, 12> mass_impl() override {
        StaticMatrix<12, 12> T = transformation();
        StaticMatrix<12, 12> M = StaticMatrix<12, 12>::Zero();

        Precision A = get_profile()->A;
        Precision L = length();
        Precision rho = get_material()->get_density();
        Precision Ip = get_profile()->I_y + get_profile()->I_z;
        Precision IpA = Ip / A;

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

        M *= rho * L * A / 420;
        return T.transpose() * M * T;
    }

    void compute_stress_strain(IPData& ip_stress,
                               IPData& ip_strain,
                               NodeData& displacement,
                               int ip_offset) override {
        (void)ip_strain;

        const Precision E = get_elasticity()->youngs;
        const Precision A = get_profile()->A;
        const Precision L = length();

        StaticMatrix<12, 1> u_global;
        for (int i = 0; i < 2; ++i) {
            Vec6 ug = displacement.row(this->nodes()[i]);
            for (int d = 0; d < 6; ++d) {
                u_global(6 * i + d) = ug(d);
            }
        }

        StaticMatrix<12, 12> T = transformation();
        StaticMatrix<12, 1> u_local = T * u_global;

        Precision eps = (u_local(6) - u_local(0)) / L;
        Precision N = E * A * eps;

        ip_stress(ip_offset, 0) = N;
        for (int j = 1; j < ip_stress.cols(); ++j) {
            ip_stress(ip_offset, j) = 0.0;
        }
    }

    StaticMatrix<12, 12> stiffness_geom_impl(IPData& ip_stress, int offset) override {
        StaticMatrix<12, 12> T = transformation();
        const Precision L = length();
        const Precision A = get_profile()->A;

        constexpr bool use_sigma = false;
        Precision Nsum = 0.0;
        int nip = 0;
        int count = 1;
        for (int i = 0; i < count; ++i) {
            const Index row = offset + i;
            Precision val = ip_stress(row, 0);
            if (use_sigma) {
                val *= A;
            }
            Nsum += val;
            ++nip;
        }
        const Precision N = (nip > 0 ? Nsum / nip : 0.0);

        if (std::abs(N) <= std::numeric_limits<Precision>::epsilon()) {
            return StaticMatrix<12, 12>::Zero();
        }

        const Precision L2 = L * L;
        const Precision f = N / (30.0 * L);

        // build Kg4 as before
        // 4×4 geometric block in local axes (already scaled by f = N/(30L))
        Eigen::Matrix<Precision, 4, 4> Kg4;
        Kg4 <<
            36.0,    3.0 * L,  -36.0,    3.0 * L,
             3.0 * L, 4.0 * L2, -3.0 * L, -1.0 * L2,
            -36.0,   -3.0 * L,  36.0,   -3.0 * L,
             3.0 * L, -1.0 * L2, -3.0 * L,  4.0 * L2;
        Kg4 *= f;

        StaticMatrix<12, 12> Kg_local = StaticMatrix<12, 12>::Zero();

        // DOF maps: [u_z, th_y] plane (“y-bending”)
        const int map_y[4] = {2, 4, 8, 10};
        // DOF maps: [u_y, th_z] plane (“z-bending”)
        const int map_z[4] = {1, 5, 7, 11};

        // Helper to scatter a 4×4 into the 12×12
        auto scatter = [&](const Eigen::Matrix<Precision,4,4>& B, const int map[4]) {
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    Kg_local(map[r], map[c]) += B(r, c);
        };

        // --- Place y-plane as-is (its (u,θ) sign already matches your K: +) ---
        scatter(Kg4, map_y);

        // --- Flip rotation sign in z-plane so (u_y, θ_z) coupling becomes negative ---
        Eigen::Matrix<Precision,4,4> S = Eigen::Matrix<Precision,4,4>::Identity();
        S(1,1) = -1; // flip θ at node 1 (second dof inside the 4×4)
        S(3,3) = -1; // flip θ at node 2 (fourth dof inside the 4×4)
        Eigen::Matrix<Precision,4,4> Kg4_z = S.transpose() * Kg4 * S;
        scatter(Kg4_z, map_z);

        return T.transpose() * Kg_local * T;
    }
};

} // namespace model
} // namespace fem

