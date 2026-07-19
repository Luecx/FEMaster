/**
 * @file c3d8r.cpp
 * @brief Implementation of the C3D8R reduced-integration solid.
 */

#include "c3d8r.h"

#include <cmath>

namespace fem::model {

// -----------------------------------------------------------------------------
// Construction and type information
// -----------------------------------------------------------------------------

C3D8R::C3D8R(ID elem_id, const std::array<ID, N>& node_ids)
    : C3D8(elem_id, node_ids) {}

std::string C3D8R::type_name() const {
    return "C3D8R";
}

// -----------------------------------------------------------------------------
// Reduced integration and output locations
// -----------------------------------------------------------------------------

const math::quadrature::Quadrature& C3D8R::integration_scheme() const {
    /*
     * One-point hexahedral rule:
     *
     *     r = s = t = 0
     *     w = 8
     */
    static const math::quadrature::Quadrature quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_CONSTANT
    };

    return quadrature;
}

RowMatrix C3D8R::stress_strain_nodal_rst() {
    return RowMatrix::Zero(N, D);
}

// -----------------------------------------------------------------------------
// Primitive hourglass modes
// -----------------------------------------------------------------------------

C3D8R::HourglassModes C3D8R::primitive_hourglass_modes() {
    /*
     * Construct the modes from the local C3D8 coordinates:
     *
     *     gamma_1 = s t
     *     gamma_2 = t r
     *     gamma_3 = r s
     *     gamma_4 = r s t.
     */
    const auto local_coords = node_coords_local();

    HourglassModes modes = HourglassModes::Zero();

    for (Index node = 0; node < N; ++node) {
        const Precision r = local_coords(node, 0);
        const Precision s = local_coords(node, 1);
        const Precision t = local_coords(node, 2);

        modes(node, 0) = s * t;
        modes(node, 1) = t * r;
        modes(node, 2) = r * s;
        modes(node, 3) = r * s * t;
    }

    return modes;
}

// -----------------------------------------------------------------------------
// Mean reference gradient
// -----------------------------------------------------------------------------

C3D8R::GradientMatrix C3D8R::mean_reference_gradient(Precision& reference_volume) {
    const auto reference_coords = node_coords_reference();

    /*
     * Do not use integration_scheme() here. The mean gradient is obtained by
     * integrating the ordinary C3D8 reference gradients with a full 2x2x2
     * rule:
     *
     *                  1
     *     D_bar = ----------- integral(D dV0).
     *                 V0
     */
    static const math::quadrature::Quadrature full_quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };

    GradientMatrix integrated_gradient = GradientMatrix::Zero();
    reference_volume                   = Precision(0);

    for (Index q = 0; q < full_quadrature.count(); ++q) {
        const auto point = full_quadrature.get_point(q);

        Precision det0 = Precision(0);
        const auto gradient = shape_derivatives_reference(
            reference_coords,
            point.r,
            point.s,
            point.t,
            det0
        );

        logging::error(std::isfinite(det0) && det0 > Precision(0),
                       "C3D8R: invalid reference determinant in element ", elem_id,
                       "\ndet(J0): ", det0);

        const Precision measure = det0 * point.w;

        integrated_gradient += gradient * measure;
        reference_volume     += measure;
    }

    logging::error(std::isfinite(reference_volume) && reference_volume > Precision(0),
                   "C3D8R: invalid reference volume in element ", elem_id,
                   "\nvolume: ", reference_volume);

    return integrated_gradient / reference_volume;
}

// -----------------------------------------------------------------------------
// Hourglass material scaling
// -----------------------------------------------------------------------------

Precision C3D8R::hourglass_material_scale() {
    /*
     * Evaluate the constitutive tangent in the undeformed reference state:
     *
     *     E = 0  <=>  F = I.
     *
     * Using the mean shear diagonal avoids coupling the hourglass stiffness
     * directly to the bulk response of nearly incompressible materials:
     *
     *     G_eff = (C44 + C55 + C66) / 3.
     */
    const Mat6 material_tangent = material_tangent_reference(
        Precision(0),
        Precision(0),
        Precision(0)
    );

    const Precision shear_scale =
        (
            material_tangent(3, 3)
            + material_tangent(4, 4)
            + material_tangent(5, 5)
        ) / Precision(3);

    logging::error(std::isfinite(shear_scale) && shear_scale > Precision(0),
                   "C3D8R: invalid initial mean shear stiffness in element ", elem_id,
                   "\nscale: ", shear_scale);

    return shear_scale;
}

// -----------------------------------------------------------------------------
// Hourglass stiffness
// -----------------------------------------------------------------------------

C3D8R::Matrix24 C3D8R::hourglass_stiffness() {
    const auto reference_coords = node_coords_reference();

    Precision reference_volume = Precision(0);

    const GradientMatrix mean_gradient = mean_reference_gradient(reference_volume);

    /*
     * Flanagan-Belytschko projection against the affine coordinate field:
     *
     *     G = (I - D_bar X^T) gamma.
     */
    const StaticMatrix<N, N> projector =
        StaticMatrix<N, N>::Identity() - mean_gradient * reference_coords.transpose();
    const HourglassModes modes = projector * primitive_hourglass_modes();

    /*
     * Constant reference scale:
     *
     *     k_hg = alpha G_eff V0 sum_a ||grad_bar N_a||^2.
     */
    const Precision material_scale  = hourglass_material_scale();
    const Precision gradient_scale  = mean_gradient.array().square().sum();
    const Precision hourglass_scale =
        default_hourglass_coefficient
        * material_scale
        * reference_volume
        * gradient_scale;

    logging::error(std::isfinite(hourglass_scale) && hourglass_scale > Precision(0),
                   "C3D8R: invalid hourglass stiffness in element ", elem_id,
                   "\nscale: ", hourglass_scale);

    /*
     * Scalar nodal matrix:
     *
     *     H = k_hg G G^T.
     */
    const StaticMatrix<N, N> scalar_stiffness =
        hourglass_scale * modes * modes.transpose();

    Matrix24 hourglass_stiffness = Matrix24::Zero();

    // Expand the scalar matrix as K_hg = kron(H, I3)
    for (Index node_a = 0; node_a < N; ++node_a) {
        for (Index node_b = 0; node_b < N; ++node_b) {
            for (Dim dof = 0; dof < D; ++dof) {
                hourglass_stiffness(D * node_a + dof, D * node_b + dof) =
                    scalar_stiffness(node_a, node_b);
            }
        }
    }

    return Precision(0.5) * (hourglass_stiffness + hourglass_stiffness.transpose());
}

// -----------------------------------------------------------------------------
// Element displacement and force assembly
// -----------------------------------------------------------------------------

C3D8R::Vector24 C3D8R::local_displacement() {
    const GradientMatrix reference_coords = node_coords_reference();
    const GradientMatrix current_coords   = node_coords_current();
    const GradientMatrix displacement     = current_coords - reference_coords;

    Vector24 result = Vector24::Zero();

    for (Index node = 0; node < N; ++node) {
        for (Dim dof = 0; dof < D; ++dof) {
            result(D * node + dof) = displacement(node, dof);
        }
    }

    return result;
}

void C3D8R::assemble_local_force(Field& node_forces, const Vector24& local_force) {
    logging::error(node_forces.domain == FieldDomain::NODE,
                   "C3D8R: internal force output must use NODE domain");
    logging::error(node_forces.components >= D,
                   "C3D8R: internal force output requires at least three components");

    for (Index node = 0; node < N; ++node) {
        const Index node_id = static_cast<Index>(node_ids[node]);

        for (Dim dof = 0; dof < D; ++dof) {
            node_forces(node_id, dof) += local_force(D * node + dof);
        }
    }
}

// -----------------------------------------------------------------------------
// Tangent stiffness
// -----------------------------------------------------------------------------

MapMatrix C3D8R::stiffness(Precision* buffer) {
    /*
     * The inherited C3D8 stiffness uses this element's virtual one-point
     * integration rule. Add the constant reference hourglass tangent:
     *
     *     K = K_cont + K_hg.
     */
    MapMatrix mapped{buffer, ndof, ndof};

    C3D8::stiffness(buffer);
    mapped += hourglass_stiffness();
    mapped  = Precision(0.5) * (mapped + mapped.transpose());

    return mapped;
}

// -----------------------------------------------------------------------------
// Nonlinear internal force
// -----------------------------------------------------------------------------

void C3D8R::compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) {
    /*
     * The inherited continuum residual uses the one-point Total-Lagrange
     * formulation with Green-Lagrange strain and PK2 stress.
     */
    C3D8::compute_internal_force_nonlinear(node_forces, ip_stress);

    /*
     * Constant reference hourglass force:
     *
     *     f_hg = K_hg u_e,
     *
     * whose exact derivative is the same constant K_hg.
     */
    const Vector24 hourglass_force = hourglass_stiffness() * local_displacement();

    assemble_local_force(node_forces, hourglass_force);
}

} // namespace fem::model
