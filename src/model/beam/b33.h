/**
 * @file b33.h
 * @brief Defines the two-node three-dimensional Euler-Bernoulli beam element.
 */

#pragma once

#include "beam.h"
#include "../geometry/line/line2a.h"

#include <limits>

namespace fem {
namespace model {

/**
 * @brief Two-node three-dimensional Euler-Bernoulli beam element.
 *
 * The element provides the linear elastic beam stiffness, consistent mass and
 * the classical initial-stress geometric stiffness used by linear buckling.
 * Prestress is derived directly from the supplied nodal displacement field;
 * no integration-point force/stress scratch field is required.
 *
 * A fully consistent finite-rotation nonlinear beam residual is not implemented
 * by B33, so nonlinear tangent evaluation remains intentionally unsupported in
 * the common `BeamElement` base.
 */
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

    void compute_stress_strain(
        Field*           strain,
        Field*           stress,
        const Field&     target_displacement,
        const Field*     target_temperature,
        const RowMatrix& rst,
        const Field*     base_displacement,
        const Field*     base_temperature
    ) override {
        // First compiled output row belonging to this element.
        Index offset = static_cast<Index>(this->elem_nodal_offset);
        if ((strain && strain->domain == FieldDomain::ELEMENT_IP) || (stress && stress->domain == FieldDomain::ELEMENT_IP)) {
            offset = static_cast<Index>(this->elem_ip_offset);
        }

        (void) base_temperature;
        logging::error(base_displacement == nullptr,
            "B33: nonlinear stress/strain evaluation is not implemented yet for element ", this->elem_id);
        logging::error(strain != nullptr || stress != nullptr,
            "B33: compute_stress_strain requires at least one output field");

        const Precision E = get_elasticity()->youngs;
        const Precision A = get_profile()->area_;
        const Precision L = length();

        StaticMatrix<12, 1> u_global;
        for (int i = 0; i < 2; ++i) {
            Vec6 ug = target_displacement.row_vec6(static_cast<Index>(this->nodes()[i]));
            for (int d = 0; d < 6; ++d) u_global(6 * i + d) = ug(d);
        }

        StaticMatrix<12, 12> T = transformation();
        StaticMatrix<12, 1> u_local = T * u_global;

        const Precision axial_strain = (u_local(6) - u_local(0)) / L;

        // Thermal expansion enters the axial constitutive response through
        //
        //     epsilon_mech = epsilon - alpha (T - T0).
        //
        // The two-node beam uses the mean nodal temperature.
        Precision thermal_strain = Precision(0);
        auto material = get_material();
        if (target_temperature && material->has_thermal_expansion()) {
            const Precision T0    = material->get_thermal_zero_temperature();
            const Precision alpha = material->get_thermal_expansion();

            Precision temperature = Precision(0);
            for (Index node = 0; node < 2; ++node) {
                temperature += (*target_temperature)(
                    static_cast<Index>(this->node_ids[node]), 0);
            }
            temperature *= Precision(0.5);
            thermal_strain = alpha * (temperature - T0);
        }

        const Precision axial_force = E * A * (axial_strain - thermal_strain);

        for (Eigen::Index i = 0; i < rst.rows(); ++i) {
            const Index row = static_cast<Index>(offset + i);
            if (strain) {
                for (Index j = 0; j < strain->components; ++j) (*strain)(row, j) = Precision(0);
                (*strain)(row, 0) = axial_strain;
            }
            if (stress) {
                for (Index j = 0; j < stress->components; ++j) (*stress)(row, j) = Precision(0);
                (*stress)(row, 0) = axial_force;
            }
        }
    }

    /**
     * Builds the classical beam-column geometric stiffness from the supplied
     * displacement state.
     *
     * The axial prestress is recovered locally as
     *
     *     Delta N = EA [(u2_x - u1_x) / L
     *                  - epsilon_th(T) + epsilon_th(T0)]
     *
     * in the beam principal frame. The resulting initial-stress matrix is then
     * rotated back to global element coordinates. No global integration-point
     * stress/resultant field is required.
     *
     * @param target_displacement Global target displacement from u0 = 0.
     * @param target_temperature Target temperature state T.
     * @param base_temperature Base temperature state T0.
     * @return Global twelve-by-twelve geometric stiffness matrix.
     */
    StaticMatrix<12, 12> stiffness_geom_impl(
        const Field& target_displacement,
        const Field* target_temperature,
        const Field* base_temperature
    ) override {
        const StaticMatrix<12, 12> T = transformation();
        const Precision L = length();

        StaticVector<12> u_global;
        for (Index node = 0; node < 2; ++node) {
            const Vec6 u = target_displacement.row_vec6(static_cast<Index>(node_ids[node]));
            for (Index dof = 0; dof < 6; ++dof) {
                u_global(6 * node + dof) = u(dof);
            }
        }

        const StaticVector<12> u_local = T * u_global;

        auto material = get_material();
        const auto thermal_strain = [&](const Field* temperature_field) {
            if (!temperature_field || !material->has_thermal_expansion()) {
                return Precision(0);
            }

            const Precision T0    = material->get_thermal_zero_temperature();
            const Precision alpha = material->get_thermal_expansion();

            Precision temperature = Precision(0);
            for (Index node = 0; node < 2; ++node) {
                temperature += (*temperature_field)(
                    static_cast<Index>(this->node_ids[node]), 0);
            }

            temperature *= Precision(0.5);
            return alpha * (temperature - T0);
        };

        const Precision N_axial =
            get_elasticity()->youngs * get_profile()->area_
            * ((u_local(6) - u_local(0)) / L
             - thermal_strain(target_temperature)
             + thermal_strain(base_temperature));

        if (std::abs(N_axial) <= std::numeric_limits<Precision>::epsilon()) {
            return StaticMatrix<12, 12>::Zero();
        }

        const Precision L2 = L * L;
        const Precision f  = N_axial / (Precision(30) * L);

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
