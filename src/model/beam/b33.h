/**
 * @file b33.h
 * @brief Defines the two-node three-dimensional Euler-Bernoulli beam element.
 */

#pragma once

#include "beam.h"
#include "../../material/strain/beam_generalized_strain.h"
#include "../geometry/line/line2a.h"

#include <limits>

namespace fem {
namespace model {

struct B33 : BeamElement<2> {
    B33(ID elem_id, std::array<ID, 2> node_ids_in)
        : BeamElement(elem_id, node_ids_in) {}

    // Instance expansion copies only persistent element topology. Sections,
    // dense ids and runtime state are assigned by Model::compile() afterwards.
    ElementPtr copy() const override { return std::make_shared<B33>(elem_id, node_ids); }

    std::string type_name() const override { return "B33"; }

    StaticMatrix<12, 12> stiffness_impl() override {
        const StaticMatrix<12, 12> Trot = transformation();
        StaticMatrix<12, 12> K_sp = StaticMatrix<12, 12>::Zero();

        Precision E   = get_elasticity()->youngs;
        Precision G   = get_elasticity()->shear;
        Precision A   = get_profile()->area_;
        Precision Iy  = get_profile()->inertia_y_;
        Precision Iz  = get_profile()->inertia_z_;
        Precision Iyz = get_profile()->product_inertia_yz_;
        Precision It  = get_profile()->torsion_inertia_;
        Precision L   = length();

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

        Profile* pr = get_profile();
        Precision ey   = pr->offset_y_;
        Precision ez   = pr->offset_z_;
        Precision refy = pr->reference_y_;
        Precision refz = pr->reference_z_;

        BeamElement<2>::rotate_yz_to_principal(phi, ey, ez);
        BeamElement<2>::rotate_yz_to_principal(phi, refy, refz);

        const StaticMatrix<12, 12> B_smp_to_sp = BeamElement<2>::rigid_offset_N(ey, ez);
        const StaticMatrix<12, 12> K_smp = B_smp_to_sp.transpose() * K_sp * B_smp_to_sp;
        const StaticMatrix<12, 12> B_ref_to_smp = BeamElement<2>::rigid_offset_N(-refy, -refz);
        const StaticMatrix<12, 12> K_ref = B_ref_to_smp.transpose() * K_smp * B_ref_to_smp;
        return Trot.transpose() * K_ref * Trot;
    }

    StaticMatrix<12, 12> mass_impl() override {
        const StaticMatrix<12, 12> Trot = transformation();
        StaticMatrix<12, 12> M_sp = StaticMatrix<12, 12>::Zero();

        Precision A   = get_profile()->area_;
        Precision L   = length();
        Precision rho = get_material()->get_density();
        Precision Ip  = get_profile()->inertia_y_ + get_profile()->inertia_z_;
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

        Precision phi = Precision(0);
        {
            const Precision Iy  = get_profile()->inertia_y_;
            const Precision Iz  = get_profile()->inertia_z_;
            const Precision Iyz = get_profile()->product_inertia_yz_;
            const Precision scale = std::max<Precision>(Precision(1), std::abs(Iy) + std::abs(Iz));
            if (std::abs(Iyz) > scale * Precision(1e-14)) phi = principal_angle();
        }

        Profile* pr = get_profile();
        Precision ey   = pr->offset_y_;
        Precision ez   = pr->offset_z_;
        Precision refy = pr->reference_y_;
        Precision refz = pr->reference_z_;

        BeamElement<2>::rotate_yz_to_principal(phi, ey, ez);
        BeamElement<2>::rotate_yz_to_principal(phi, refy, refz);

        const StaticMatrix<12, 12> B_smp_to_sp = BeamElement<2>::rigid_offset_N(ey, ez);
        const StaticMatrix<12, 12> M_smp = B_smp_to_sp.transpose() * M_sp * B_smp_to_sp;
        const StaticMatrix<12, 12> B_ref_to_smp = BeamElement<2>::rigid_offset_N(-refy, -refz);
        const StaticMatrix<12, 12> M_ref = B_ref_to_smp.transpose() * M_smp * B_ref_to_smp;
        return Trot.transpose() * M_ref * Trot;
    }

    RowMatrix stress_strain_ip_rst() override {
        RowMatrix rst(1, 3);
        rst.setZero();
        return rst;
    }

    void compute_stress_strain(Field* strain,
                               Field* stress,
                               const Field& displacement,
                               const RowMatrix& rst,
                               int offset,
                               bool use_green_lagrange_nl) override {
        logging::error(!use_green_lagrange_nl,
            "B33: nonlinear stress/strain evaluation is not implemented yet for element ", this->elem_id);
        logging::error(strain != nullptr || stress != nullptr,
            "B33: compute_stress_strain requires at least one output field");

        const Precision E = get_elasticity()->youngs;
        const Precision A = get_profile()->area_;
        const Precision L = length();

        StaticMatrix<12, 1> u_global;
        for (int i = 0; i < 2; ++i) {
            Vec6 ug = displacement.row_vec6(static_cast<Index>(this->nodes()[i]));
            for (int d = 0; d < 6; ++d) u_global(6 * i + d) = ug(d);
        }

        StaticMatrix<12, 12> T = transformation();
        StaticMatrix<12, 1> u_local = T * u_global;

        Vec6 generalized_values = Vec6::Zero();
        generalized_values(0) = (u_local(6) - u_local(0)) / L;
        const BeamGeneralizedStrain generalized_strain(generalized_values);
        const Precision axial_force = E * A * generalized_strain.values()(0);

        for (Eigen::Index i = 0; i < rst.rows(); ++i) {
            const Index row = static_cast<Index>(offset + i);
            if (strain) {
                for (Index j = 0; j < strain->components; ++j) (*strain)(row, j) = Precision(0);
                (*strain)(row, 0) = generalized_strain.values()(0);
            }
            if (stress) {
                for (Index j = 0; j < stress->components; ++j) (*stress)(row, j) = Precision(0);
                (*stress)(row, 0) = axial_force;
            }
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
              3.0 * L,  4.0 * L2, -3.0 * L, -1.0 * L2,
            -36.0    , -3.0 * L ,  36.0    , -3.0 * L,
              3.0 * L, -1.0 * L2, -3.0 * L,  4.0 * L2;
        Kg41 *= f;
        Kg42 <<
             36.0    , -3.0 * L , -36.0    , -3.0 * L,
             -3.0 * L,  4.0 * L2,   3.0 * L, -1.0 * L2,
            -36.0    ,  3.0 * L ,  36.0    ,  3.0 * L,
             -3.0 * L, -1.0 * L2,   3.0 * L,  4.0 * L2;
        Kg42 *= f;

        StaticMatrix<12, 12> Kg_local = StaticMatrix<12, 12>::Zero();
        const int map_y[4] = {2, 4, 8, 10};
        const int map_z[4] = {1, 5, 7, 11};

        auto scatter = [&](const Eigen::Matrix<Precision, 4, 4>& B, const int map[4]) {
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    Kg_local(map[r], map[c]) += B(r, c);
        };

        scatter(Kg42, map_y);
        scatter(Kg41, map_z);
        return T.transpose() * Kg_local * T;
    }

    LinePtr line(ID line_id) override {
        (void) line_id;
        return std::make_shared<Line2A>(this->node_ids);
    }
};

} // namespace model
} // namespace fem
