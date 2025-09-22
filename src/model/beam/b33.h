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

    // Store axial force N into ip_stress at given offset.
    // We assume ip_stress has 6 columns like solids: [σxx, σyy, σzz, σyz, σzx, σxy].
    // For beams, we only fill column 0 with axial force resultant N.
    void compute_stress_strain(IPData& ip_stress,
                               IPData& ip_strain,
                               NodeData& displacement,
                               int ip_offset) override {
        (void) ip_strain; // not used for beams here

        // ---- Geometry and material data ----
        const Precision E  = get_elasticity()->youngs;
        const Precision A  = get_profile()->A;
        const Precision L  = length();

        // ---- Get nodal axial displacements in local coords ----
        // Transform global displacements to local system
        StaticMatrix<12,1> u_global;
        for (int i = 0; i < 2; ++i) {
            Vec6 ug = displacement.row(this->nodes()[i]);
            for (int d = 0; d < 6; ++d)
                u_global(6*i + d) = ug(d);
        }

        // Apply transformation to local system
        StaticMatrix<12,12> T = transformation();
        StaticMatrix<12,1> u_local = T * u_global;

        // ---- Axial force (resultant) ----
        // Axial strain = (u2x - u1x)/L
        Precision eps = (u_local(6) - u_local(0)) / L;
        Precision N   = E * A * eps; // axial force resultant

        // ---- Store into ip_stress ----
        // Here we only fill σxx column with N (resultant).
        // If you prefer to store stress (σ = N/A), replace with (N / A).
        ip_stress(ip_offset, 0) = N;

        // Other entries left zero
        for (int j = 1; j < ip_stress.cols(); ++j) {
            ip_stress(ip_offset, j) = 0.0;
        }
    }


    StaticMatrix<12,12> stiffness_geom_impl(IPData& ip_stress, int offset) override {
        // Local transform and geometry
        StaticMatrix<12,12> T = transformation();
        const Precision L = length();
        const Precision A = get_profile()->A;

        // ---- 1) Get axial force N from ip_stress (average over this element’s IPs) ----
        // Assumes ip_stress(rows, 0) = either axial resultant N or axial stress sigma_xx.
        // If you store sigma_xx: set 'use_sigma' = true to multiply by A.
        constexpr bool use_sigma = false; // set true if ip_stress contains sigma_xx (not N)
        Precision Nsum = 0.0;
        int nip = 0;
        // If your beam has its own scheme(), use it. Otherwise, assume contiguous rows per element.
        // Here we assume your integration scheme count equals number of rows you wrote for this element:
        // (If you know your exact count, use it; otherwise, you can fall back to 1 IP.)
        int count = 1; // default if you only store one IP; change if you store more
        // If you DID store multiple IPs, replace the next two lines by your real 'count'.
        // e.g., int count = this->integration_scheme().count();

        for (int i = 0; i < count; ++i) {
            const Index row = offset + i;
            Precision val = ip_stress(row, 0); // either N or sigma_xx
            if (use_sigma) val *= A;           // convert stress -> resultant
            Nsum += val;
            ++nip;
        }
        const Precision N = (nip > 0 ? Nsum / nip : 0.0); // axial force resultant (tension +, compression −)

        // Early out: no axial force -> no geometric stiffness
        if (std::abs(N) <= std::numeric_limits<Precision>::epsilon()) {
            return StaticMatrix<12,12>::Zero();
        }

        // ---- 2) Build 4x4 geometric submatrix Kg4(N,L) ----
        const Precision L2 = L*L;
        const Precision f  = N / (30.0 * L); // scalar factor

        Eigen::Matrix<Precision,4,4> Kg4;
        Kg4 <<  36.0,   3.0*L,  -36.0,   3.0*L,
                 3.0*L,  4.0*L2, -3.0*L, -1.0*L2,
                -36.0,  -3.0*L,  36.0,  -3.0*L,
                 3.0*L, -1.0*L2, -3.0*L,  4.0*L2;
        Kg4 *= f;

        // ---- 3) Assemble into 12x12 local matrix ----
        StaticMatrix<12,12> Kg_local = StaticMatrix<12,12>::Zero();

        auto add_block = [&](int i0, int j0, const Eigen::Matrix<Precision,4,4>& B){
            // i0/j0 are DOF indices [0..11] for the ordered set [v1, r1, v2, r2]
            const int map[4] = { i0, j0, i0+6, j0+6 }; // careful: not simply +6 for rotations; we pass exact indices below
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    Kg_local(map[r], map[c]) += B(r,c);
        };

        // Bending in the x–y plane (deflection uz, rotation ry): DOFs [uz1, ry1, uz2, ry2] -> indices [2,4,8,10]
        {
            const int I[4] = { 2, 4, 8, 10 };
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    Kg_local(I[r], I[c]) += Kg4(r,c);
        }

        // Bending in the x–z plane (deflection uy, rotation rz): DOFs [uy1, rz1, uy2, rz2] -> indices [1,5,7,11]
        {
            const int I[4] = { 1, 5, 7, 11 };
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    Kg_local(I[r], I[c]) += Kg4(r,c);
        }

        // (No contribution to axial ux or torsion rx in the classical beam-column Kg.)

        // ---- 4) Transform back to global ----
        return T.transpose() * Kg_local * T;
    }

};

}