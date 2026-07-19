/**
 * @file b33.h
 * @brief Defines the two-node three-dimensional Euler-Bernoulli beam element.
 *
 * `B33` specializes `BeamElement<2>` for a straight two-node beam with six
 * degrees of freedom per node. The element provides the elastic stiffness,
 * consistent mass and initial-stress geometric stiffness matrices, including
 * transformations between the global frame, the principal section frame and
 * offset reference axes.
 *
 * The current stress recovery exposes the constant axial generalized strain
 * and axial force at one integration point. Nonlinear stress recovery and a
 * nonlinear internal-force formulation are not implemented by this element.
 *
 * @see BeamElement
 * @see beam.h
 */

#pragma once

#include "beam.h"
#include "../../material/strain/beam_generalized_strain.h"
#include "../geometry/line/line2a.h"

#include <limits>

namespace fem {
namespace model {

/**
 * @brief Two-node spatial Euler-Bernoulli beam element with section offsets.
 *
 * The element uses twelve nodal degrees of freedom ordered as three
 * translations and three rotations at each node. Elastic bending is assembled
 * in the principal section frame, while rigid-offset transformations transfer
 * stiffness and inertia from the section-property axis through the section
 * midpoint to the element reference line before the final global rotation.
 *
 * The geometric stiffness is assembled from the recovered axial force, and the
 * line representation exposes the two-node element centerline for generic
 * geometry operations.
 */
struct B33 : BeamElement<2> {
    // Construction from the element identifier, the two ordered beam nodes and
    // an optional orientation node. All common state is initialized by the
    // `BeamElement<2>` base class.
    B33(ID elem_id, std::array<ID, 2> node_ids_in, ID orientation_node_id = static_cast<ID>(-1))
        : BeamElement(elem_id, node_ids_in, orientation_node_id) {}

    // Return the element keyword used by model input, diagnostics and output.
    std::string type_name() const override { return "B33"; }

    // Assemble the global elastic stiffness matrix. The classical uncoupled
    // Euler-Bernoulli matrix is formed about the section-property axis in the
    // principal bending frame, transferred through both section offsets and
    // finally rotated into global degrees of freedom.
    StaticMatrix<12, 12> stiffness_impl() override {
        // Transform global nodal degrees of freedom to the principal section frame.
        const StaticMatrix<12, 12> Trot = transformation();

        // ------------------------------------------------------------------
        // 1) Form the uncoupled Euler-Bernoulli stiffness about the
        //    section-property axis in the principal bending frame.
        // ------------------------------------------------------------------
        StaticMatrix<12, 12> K_sp = StaticMatrix<12, 12>::Zero();

        Precision E   = get_elasticity()->youngs;
        Precision G   = get_elasticity()->shear;
        Precision A   = get_profile()->area_;
        Precision Iy  = get_profile()->inertia_y_;
        Precision Iz  = get_profile()->inertia_z_;
        Precision Iyz = get_profile()->product_inertia_yz_;
        Precision It  = get_profile()->torsion_inertia_;
        Precision L   = length();

        // Diagonalize the bending inertia tensor when a relevant product
        // inertia is present. Section offsets are rotated by the same angle below.
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
            // In the rotated basis the product inertia is zero up to the
            // scale-aware numerical tolerance used for the principal-axis test.
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
        // 2) Read the offsets from the section midpoint to the section-property
        //    axis and to the element reference axis. Because K_sp is expressed
        //    in principal coordinates, both offsets are rotated into the same
        //    local y-z frame before the kinematic transfers are applied.
        // ------------------------------------------------------------------
        Profile* pr = get_profile();

        // Offsets are stored relative to SMP and are rotated into the same
        // local basis as the principal section stiffness before assembly.
        Precision ey   = pr->offset_y_;
        Precision ez   = pr->offset_z_;
        Precision refy = pr->reference_y_;
        Precision refz = pr->reference_z_;

        // Rotate the stored base-frame offsets into principal coordinates.
        BeamElement<2>::rotate_yz_to_principal(phi, ey,   ez);
        BeamElement<2>::rotate_yz_to_principal(phi, refy, refz);

        // ------------------------------------------------------------------
        // 3) Transfer the section-property axis stiffness to the section midpoint.
        //
        //    u_SP = B(SP - SMP) * u_SMP
        //    K_SMP = B^T * K_SP * B
        // ------------------------------------------------------------------
        const StaticMatrix<12, 12> B_smp_to_sp = BeamElement<2>::rigid_offset_N(ey, ez);
        const StaticMatrix<12, 12> K_smp = B_smp_to_sp.transpose() * K_sp * B_smp_to_sp;

        // ------------------------------------------------------------------
        // 4) Transfer the section-midpoint stiffness to the reference axis.
        //
        //    ref_y,z = REF - SMP  =>  SMP - REF = -(REF - SMP) = (-refy, -refz)
        //
        //    u_SMP = B(SMP - REF) * u_REF
        //    K_REF = B^T * K_SMP * B
        // ------------------------------------------------------------------
        const StaticMatrix<12, 12> B_ref_to_smp = BeamElement<2>::rigid_offset_N(-refy, -refz);
        const StaticMatrix<12, 12> K_ref = B_ref_to_smp.transpose() * K_smp * B_ref_to_smp;

        // ------------------------------------------------------------------
        // 5) Rotate the reference-axis matrix from principal local to global DOFs.
        // ------------------------------------------------------------------
        return Trot.transpose() * K_ref * Trot;
    }

    // Assemble the global consistent mass matrix. The local matrix includes
    // translational inertia, rotary inertia about the beam axis and the same
    // two-stage rigid-offset mapping used by the elastic stiffness.
    StaticMatrix<12, 12> mass_impl() override {
        // Transform global nodal degrees of freedom to the principal section frame.
        const StaticMatrix<12, 12> Trot = transformation();

        // ------------------------------------------------------------------
        // 1) Form the classical consistent mass matrix about the
        //    section-property axis.
        // ------------------------------------------------------------------
        StaticMatrix<12, 12> M_sp = StaticMatrix<12, 12>::Zero();

        Precision A   = get_profile()->area_;
        Precision L   = length();
        Precision rho = get_material()->get_density();
        // The rotary inertia term uses the polar inertia about the
        // section-property axis, obtained as the sum of the two transverse
        // second moments.
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

        // ------------------------------------------------------------------
        // 2) Reuse the principal-axis angle employed by the degree-of-freedom
        //    transformation and express all offsets in that frame.
        // ------------------------------------------------------------------
        Precision phi = Precision(0);
        {
            const Precision Iy  = get_profile()->inertia_y_;
            const Precision Iz  = get_profile()->inertia_z_;
            const Precision Iyz = get_profile()->product_inertia_yz_;
            const Precision scale = std::max<Precision>(Precision(1), std::abs(Iy) + std::abs(Iz));
            if (std::abs(Iyz) > scale * Precision(1e-14)) {
                phi = principal_angle();
            }
        }

        Profile* pr = get_profile();

        // Offsets from the section midpoint to the section-property axis
        // and to the element reference axis, expressed initially in the base frame.
        Precision ey   = pr->offset_y_;
        Precision ez   = pr->offset_z_;
        Precision refy = pr->reference_y_;
        Precision refz = pr->reference_z_;

        BeamElement<2>::rotate_yz_to_principal(phi, ey,   ez);
        BeamElement<2>::rotate_yz_to_principal(phi, refy, refz);

        // ------------------------------------------------------------------
        // 3) Apply the same section-property-axis to midpoint and midpoint
        //    to reference-axis mappings used by the stiffness matrix.
        // ------------------------------------------------------------------
        const StaticMatrix<12, 12> B_smp_to_sp  = BeamElement<2>::rigid_offset_N(ey, ez);
        const StaticMatrix<12, 12> M_smp        = B_smp_to_sp.transpose() * M_sp * B_smp_to_sp;

        // Map the reference-axis DOFs to the section-midpoint DOFs using the
        // negative of the stored REF-SMP offset.
        const StaticMatrix<12, 12> B_ref_to_smp = BeamElement<2>::rigid_offset_N(-refy, -refz);
        const StaticMatrix<12, 12> M_ref        = B_ref_to_smp.transpose() * M_smp * B_ref_to_smp;

        // ------------------------------------------------------------------
        // 4) Rotate the reference-axis mass matrix into global coordinates.
        // ------------------------------------------------------------------
        return Trot.transpose() * M_ref * Trot;
	}

    // Return the single beam integration point used by the current constant
    // axial stress and strain recovery. Its natural coordinates are all zero.
    RowMatrix stress_strain_ip_rst() override {
        RowMatrix rst(1, 3);
        rst.setZero();
        return rst;
    }

    // Recover the constant axial generalized strain and corresponding axial
    // force from the local end displacement difference. Requested result fields
    // are cleared at every supplied integration-point row before component zero
    // is populated.
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
            for (int d = 0; d < 6; ++d) {
                u_global(6 * i + d) = ug(d);
            }
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

    // Assemble the global initial-stress geometric stiffness from the axial
    // force stored at the selected integration-point row. The two bending-plane
    // submatrices are scattered into the twelve local beam degrees of freedom
    // and then transformed to global coordinates.
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

        // Scatter the first bending-plane matrix into z-translation/theta-y
        // DOFs and the second into y-translation/theta-z DOFs.
        const int map_y[4] = {2, 4, 8, 10};
        const int map_z[4] = {1, 5, 7, 11};

        auto scatter = [&](const Eigen::Matrix<Precision,4,4>& B, const int map[4]) {
            for (int r = 0; r < 4; ++r)
                for (int c = 0; c < 4; ++c)
                    Kg_local(map[r], map[c]) += B(r, c);
        };

        scatter(Kg42, map_y);
        scatter(Kg41, map_z);

        return T.transpose() * Kg_local * T;
    }

    // Create the two-node geometric centerline representation. A B33 element
    // exposes only one line, so the requested local line identifier is ignored.
    LinePtr line(ID line_id) override {
        (void)line_id;
        return std::make_shared<Line2A>(this->node_ids);
    }
};
} // namespace model
} // namespace fem
