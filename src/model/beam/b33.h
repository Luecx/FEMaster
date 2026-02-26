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
#include "../geometry/line/line2a.h"

#include <limits>

namespace fem {
namespace model {

struct B33 : BeamElement<2> {
    B33(ID elem_id, std::array<ID, 2> node_ids_in, ID orientation_node_id = static_cast<ID>(-1))
        : BeamElement(elem_id, node_ids_in, orientation_node_id) {}

    std::string type_name() const override { return "B33"; }

    StaticMatrix<12, 12> stiffness_impl() override {
        // This T maps global DOFs -> local DOFs in the PRINCIPAL frame (see BeamElement::transformation()).
        const StaticMatrix<12, 12> Trot = transformation();

        // ------------------------------------------------------------------
        // 1) Build classical Euler-Bernoulli K_loc referenced to SP-axis
        //    (i.e. the standard uncoupled beam stiffness).
        // ------------------------------------------------------------------
        StaticMatrix<12, 12> K_sp = StaticMatrix<12, 12>::Zero();

        Precision E = get_elasticity()->youngs;
        Precision G = get_elasticity()->shear;
        Precision A = get_profile()->A;
        Precision Iy = get_profile()->I_y;
        Precision Iz = get_profile()->I_z;
        Precision Iyz = get_profile()->I_yz;
        Precision It = get_profile()->I_t;
        Precision L = length();

        // rotate inertias to principal axes (and we will rotate offsets with the same phi)
        Precision phi = Precision(0);
        const Precision scale = std::max<Precision>(Precision(1), std::abs(Iy) + std::abs(Iz));
        if (std::abs(Iyz) > scale * Precision(1e-14)) {
            phi = principal_angle();
            const Precision cph = std::cos(phi);
            const Precision sph = std::sin(phi);
            const Precision c2 = cph * cph;
            const Precision s2 = sph * sph;
            const Precision sc = sph * cph;
            const Precision Iy_p = Iy * c2 + Iz * s2 - 2 * Iyz * sc;
            const Precision Iz_p = Iy * s2 + Iz * c2 + 2 * Iyz * sc;
            Iy = Iy_p;
            Iz = Iz_p;
            // (Iyz becomes ~0 in this basis)
        }

        const Precision a = E * A / L;
        const Precision b = G * It / L;
        const Precision c = E * Iz / (L * L * L);
        const Precision d = E * Iy / (L * L * L);
        const Precision M = L * L;

        K_sp <<
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

        // ------------------------------------------------------------------
        // 2) Offsets provided relative to SMP:
        //    ey, ez   := COG - SMP = SP - SMP
        //    ref_y,z  := REF - SMP
        //
        //    We must express these in the same LOCAL frame as K_sp:
        //    since Trot uses principal axes, rotate offsets by phi.
        // ------------------------------------------------------------------
        Profile* pr = get_profile();

        // REQUIRED profile members (rename if needed):
        //   pr->e_y, pr->e_z     : (SP - SMP)
        //   pr->ref_y, pr->ref_z : (REF - SMP)
        Precision ey = pr->e_y;
        Precision ez = pr->e_z;
        Precision refy = pr->ref_y;
        Precision refz = pr->ref_z;

        // rotate offsets from base yz to principal yz if needed
        BeamElement<2>::rotate_yz_to_principal(phi, ey,   ez);
        BeamElement<2>::rotate_yz_to_principal(phi, refy, refz);

        // ------------------------------------------------------------------
        // 3) Stage 1: SP <- SMP via kinematics
        //
        //    u_SP = B(SP - SMP) * u_SMP
        //    K_SMP = B^T * K_SP * B
        // ------------------------------------------------------------------
        const StaticMatrix<12, 12> B_smp_to_sp = BeamElement<2>::rigid_offset_N(ey, ez);
        const StaticMatrix<12, 12> K_smp = B_smp_to_sp.transpose() * K_sp * B_smp_to_sp;

        // ------------------------------------------------------------------
        // 4) Stage 2: SMP <- REF via kinematics
        //
        //    ref_y,z = REF - SMP  =>  SMP - REF = -(REF - SMP) = (-refy, -refz)
        //
        //    u_SMP = B(SMP - REF) * u_REF
        //    K_REF = B^T * K_SMP * B
        // ------------------------------------------------------------------
        const StaticMatrix<12, 12> B_ref_to_smp = BeamElement<2>::rigid_offset_N(-refy, -refz);
        const StaticMatrix<12, 12> K_ref = B_ref_to_smp.transpose() * K_smp * B_ref_to_smp;

        // ------------------------------------------------------------------
        // 5) Rotate local(ref, principal) -> global
        // ------------------------------------------------------------------
        return Trot.transpose() * K_ref * Trot;
    }

    StaticMatrix<12, 12> mass_impl() override {
        // DOF rotation: global -> local (principal frame)
        const StaticMatrix<12, 12> Trot = transformation();

        // ------------------------------------------------------------------
        // 1) Build your current consistent mass matrix M_sp (as-is).
        //    (Interpreted as referenced to SP-axis, consistent with classic theory.)
        // ------------------------------------------------------------------
        StaticMatrix<12, 12> M_sp = StaticMatrix<12, 12>::Zero();

        Precision A   = get_profile()->A;
        Precision L   = length();
        Precision rho = get_material()->get_density();
        Precision Ip  = get_profile()->I_y + get_profile()->I_z; // polar about SP if I_y, I_z about SP
        Precision IpA = Ip / A;

        M_sp <<
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

        M_sp *= rho * L * A / 420;

        // ------------------------------------------------------------------
        // 2) Compute the same principal-axis rotation angle used by transformation(),
        //    and rotate offsets into the principal yz-frame.
        // ------------------------------------------------------------------
        Precision phi = Precision(0);
        {
            const Precision Iy  = get_profile()->I_y;
            const Precision Iz  = get_profile()->I_z;
            const Precision Iyz = get_profile()->I_yz;
            const Precision scale = std::max<Precision>(Precision(1), std::abs(Iy) + std::abs(Iz));
            if (std::abs(Iyz) > scale * Precision(1e-14)) {
                phi = principal_angle();
            }
        }

        Profile* pr = get_profile();

        // Offsets provided relative to SMP:
        //   ey, ez   := COG - SMP = SP - SMP
        //   ref_y,z  := REF - SMP
        Precision ey   = pr->e_y;
        Precision ez   = pr->e_z;
        Precision refy = pr->ref_y;
        Precision refz = pr->ref_z;

        BeamElement<2>::rotate_yz_to_principal(phi, ey,   ez);
        BeamElement<2>::rotate_yz_to_principal(phi, refy, refz);

        // ------------------------------------------------------------------
        // 3) Apply the exact same two-stage kinematic mapping:
        //    SP <- SMP  and  SMP <- REF
        // ------------------------------------------------------------------
        const StaticMatrix<12, 12> B_smp_to_sp  = BeamElement<2>::rigid_offset_N(ey, ez);
        const StaticMatrix<12, 12> M_smp        = B_smp_to_sp.transpose() * M_sp * B_smp_to_sp;

        const StaticMatrix<12, 12> B_ref_to_smp = BeamElement<2>::rigid_offset_N(-refy, -refz); // SMP-REF
        const StaticMatrix<12, 12> M_ref        = B_ref_to_smp.transpose() * M_smp * B_ref_to_smp;

        // ------------------------------------------------------------------
        // 4) Rotate to global
        // ------------------------------------------------------------------
        return Trot.transpose() * M_ref * Trot;
	}

    void compute_stress_strain(Field& ip_stress,
                               Field& ip_strain,
                               Field& displacement,
                               int ip_offset) override {
        (void)ip_strain;

        const Precision E = get_elasticity()->youngs;
        const Precision A = get_profile()->A;
        const Precision L = length();

        StaticMatrix<12, 1> u_global;
        for (int i = 0; i < 2; ++i) {
            Vec6 ug = displacement.row_vec6(static_cast<Index>(this->nodes()[i]));
            for (int d = 0; d < 6; ++d) {
                u_global(6 * i + d) = ug(d);
            }
        }

        // u_local is in principal local frame at REF line (because stiffness_impl maps DOFs at REF line)
        StaticMatrix<12, 12> T = transformation();
        StaticMatrix<12, 1> u_local = T * u_global;

        Precision eps = (u_local(6) - u_local(0)) / L;
        Precision N = E * A * eps;

        const Index row = static_cast<Index>(ip_offset);
        ip_stress(row, 0) = N;
        for (int j = 1; j < ip_stress.components; ++j) {
            ip_stress(row, j) = 0.0;
        }
    }

    StaticMatrix<12, 12> stiffness_geom_impl(const Field& ip_stress, int offset) override {
        StaticMatrix<12, 12> T = transformation();
        const Precision L = length();

        Precision N = ip_stress(static_cast<Index>(offset), 0);

        if (std::abs(N) <= std::numeric_limits<Precision>::epsilon()) {
            return StaticMatrix<12, 12>::Zero();
        }

        const Precision L2 = L * L;
        const Precision f = N / (30.0 * L);

        Eigen::Matrix<Precision, 4, 4> Kg41;
        Eigen::Matrix<Precision, 4, 4> Kg42;
        Kg41 <<
             36.0    ,  3.0 * L , -36.0    ,  3.0 * L,
              3.0 * L,  4.0 * L2, - 3.0 * L, -1.0 * L2,
            -36.0    , -3.0 * L ,  36.0    , -3.0 * L,
              3.0 * L, -1.0 * L2, - 3.0 * L,  4.0 * L2;
        Kg41 *= f;
        Kg42 <<
             36.0    , -3.0 * L , -36.0    , -3.0 * L,
            - 3.0 * L,  4.0 * L2,   3.0 * L, -1.0 * L2,
            -36.0    ,  3.0 * L ,  36.0    ,  3.0 * L,
            - 3.0 * L, -1.0 * L2,   3.0 * L,  4.0 * L2;
        Kg42 *= f;

        StaticMatrix<12, 12> Kg_local = StaticMatrix<12, 12>::Zero();

        const int map_y[4] = {2, 4, 8, 10};  // [u_z, th_y]
        const int map_z[4] = {1, 5, 7, 11};  // [u_y, th_z]

        auto scatter = [&](const Eigen::Matrix<Precision,4,4>& B, const int map[4]) {
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    Kg_local(map[r], map[c]) += B(r, c);
        };

        scatter(Kg42, map_y);
        scatter(Kg41, map_z);

        return T.transpose() * Kg_local * T;
    }

    LinePtr line(ID line_id) override {
        (void)line_id;
        return std::make_shared<Line2A>(this->node_ids);
    }
};

} // namespace model
} // namespace fem