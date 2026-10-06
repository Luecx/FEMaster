/**
 * @file c3d8r.cpp
 * @brief Implements the reduced-integration C3D8 solid and hourglass stabilization.
 *
 * The one-point continuum contribution uses the common solid material-point
 * state. Hourglass stiffness is obtained from a zero-strain constitutive tangent
 * without changing that physical history state, and the same hourglass matrix
 * supplies the residual and tangent contributions.
 *
 * Natural-space polynomial recovery is cached in extrapolation_matrix();
 * the common solid formulation applies it to integration-point output.
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "c3d8r.h"

#include "../../math/extrapolate.h"

#include <cmath>

namespace fem::model {

C3D8R::C3D8R(ID elem_id, const std::array<ID, N>& node_ids)
    : C3D8(elem_id, node_ids) {}

std::string C3D8R::type_name() const {
    return "C3D8R";
}

/**
 * Returns the static quadrature used by constitutive integration.
 *
 * The point ordering defines material-state rows and integration-point output.
 * This rule may contain fewer points than the topology volume quadrature.
 *
 * @return Shared immutable material quadrature.
 */
const math::quadrature::Quadrature& C3D8R::integration_scheme_stiffness() const {
    // One-point hexahedral integration at r = s = t = 0 with weight eight.
    static const math::quadrature::Quadrature quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_CONSTANT
    };

    return quadrature;
}

RowMatrix C3D8R::stress_strain_nodal_rst() {
    return RowMatrix::Zero(N, D);
}

/**
 * Constructs the four primitive scalar hourglass modes at the C3D8 nodes.
 *
 * The natural-coordinate products `[s t, t r, r s, r s t]` span the non-affine
 * zero-energy patterns of one-point hexahedral integration before projection
 * against the physical affine coordinate field.
 *
 * @return Eight-by-four primitive modal matrix in element-node ordering.
 */
C3D8R::HourglassModes C3D8R::primitive_hourglass_modes() {
    // Primitive scalar modes gamma = [st, tr, rs, rst].
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

/**
 * Computes the volume-averaged reference shape-function gradients.
 *
 * A full two-by-two-by-two rule integrates the ordinary C3D8 gradients over the
 * undeformed element:
 *
 * D_bar = (1 / V0) integral_A0 D dV0.
 *
 * Positive finite point determinants and total reference volume are required.
 * The result depends only on reference geometry and is independent of the
 * one-point continuum quadrature used by the reduced element.
 *
 * @param reference_volume Physical reference volume returned to the caller.
 * @return Mean reference gradient with one row per element node.
 */
C3D8R::GradientMatrix C3D8R::mean_reference_gradient(Precision& reference_volume) {
    const auto reference_coords = node_coords_reference();

    // Integrate the ordinary C3D8 reference gradients with a full 2x2x2 rule:
    //
    //     D_bar = (1 / V0) integral(D dV0).
    static const math::quadrature::Quadrature full_quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };

    GradientMatrix integrated_gradient = GradientMatrix::Zero();
    reference_volume                   = Precision(0);

    for (Index q = 0; q < full_quadrature.count(); ++q) {
        const auto point = full_quadrature.get_point(q);

        Precision det0      = Precision(0);
        const auto gradient = shape_derivatives_reference(reference_coords, point.r, point.s, point.t, det0);

        logging::error(std::isfinite(det0) && det0 > Precision(0),
            "C3D8R: invalid reference determinant in element ", elem_id, "\ndet(J0): ", det0);

        const Precision measure = det0 * point.w;
        integrated_gradient     += gradient * measure;
        reference_volume        += measure;
    }

    logging::error(std::isfinite(reference_volume) && reference_volume > Precision(0),
        "C3D8R: invalid reference volume in element ", elem_id, "\nvolume: ", reference_volume);

    return integrated_gradient / reference_volume;
}

/**
 * Evaluates the zero-strain constitutive shear scale without advancing the
 * physical material-point history.
 *
 * The hourglass modulus is an auxiliary stabilization quantity rather than a
 * constitutive update of the current continuum state. The zero-strain tangent
 * therefore reads the immutable committed row and supplies no target state.
 *
 * @return Mean material shear diagonal `(C44 + C55 + C66) / 3`.
 */
Precision C3D8R::hourglass_material_scale() {
    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);

    const Mat6 material_tangent = material_tangent_reference(
        Precision(0), Precision(0), Precision(0), old_state, nullptr);

    const Precision shear_scale =
        (material_tangent(3, 3) + material_tangent(4, 4) + material_tangent(5, 5)) / Precision(3);

    logging::error(std::isfinite(shear_scale) && shear_scale > Precision(0),
        "C3D8R: invalid initial mean shear stiffness in element ", elem_id, "\nscale: ", shear_scale);

    return shear_scale;
}

/**
 * Builds the constant reference hourglass stabilization tangent.
 *
 * Primitive modes are projected with `I - D_bar X^T` so affine displacement
 * fields remain unstabilized. Their scalar stiffness uses the initial mean
 * constitutive shear diagonal, reference volume and mean-gradient norm:
 *
 * k_hg = alpha G_eff V0 sum_a ||grad_bar N_a||^2.
 *
 * The scalar nodal matrix `k_hg G G^T` is expanded independently into all three
 * translational directions. The final matrix is symmetrized only to remove
 * round-off asymmetry.
 *
 * @return Constant 24-by-24 hourglass tangent in node-major DOF ordering.
 */
C3D8R::Matrix24 C3D8R::hourglass_stiffness() {
    const auto reference_coords = node_coords_reference();

    Precision reference_volume         = Precision(0);
    const GradientMatrix mean_gradient = mean_reference_gradient(reference_volume);

    // Flanagan-Belytschko projection against the affine coordinate field:
    //
    //     G = (I - D_bar X^T) gamma.
    const StaticMatrix<N, N> projector = StaticMatrix<N, N>::Identity() - mean_gradient * reference_coords.transpose();
    const HourglassModes modes         = projector * primitive_hourglass_modes();

    // Reference stabilization scale
    //
    //     k_hg = alpha G_eff V0 sum_a ||grad_bar N_a||^2.
    const Precision material_scale = hourglass_material_scale();
    const Precision gradient_scale = mean_gradient.array().square().sum();
    const Precision hourglass_scale =
        default_hourglass_coefficient * material_scale * reference_volume * gradient_scale;

    logging::error(std::isfinite(hourglass_scale) && hourglass_scale > Precision(0),
        "C3D8R: invalid hourglass stiffness in element ", elem_id, "\nscale: ", hourglass_scale);

    const StaticMatrix<N, N> scalar_stiffness = hourglass_scale * modes * modes.transpose();

    Matrix24 stiffness = Matrix24::Zero();

    // Expand the scalar matrix as K_hg = kron(H, I3).
    for (Index node_a = 0; node_a < N; ++node_a) {
        for (Index node_b = 0; node_b < N; ++node_b) {
            for (Dim dof = 0; dof < D; ++dof) {
                stiffness(D * node_a + dof, D * node_b + dof) = scalar_stiffness(node_a, node_b);
            }
        }
    }

    return Precision(0.5) * (stiffness + stiffness.transpose());
}

/**
 * Collects the supplied element translations in node-major XYZ ordering.
 *
 * The nonlinear structural API passes the trial displacement explicitly, so the
 * hourglass residual uses the same supplied state as the continuum evaluation
 * rather than reconstructing displacement from persistent current positions.
 *
 * @param displacement Global nodal trial displacement field.
 * @return Twenty-four-component displacement vector in node-major XYZ ordering.
 */
C3D8R::Vector24 C3D8R::local_displacement(const Field& displacement) {
    const GradientMatrix local = this->nodal_data<D>(displacement);
    Vector24 result            = Vector24::Zero();

    for (Index node = 0; node < N; ++node) {
        for (Dim dof = 0; dof < D; ++dof) {
            result(D * node + dof) = local(node, dof);
        }
    }

    return result;
}

/**
 * Scatters one element-local translational force vector into the global nodal
 * force field.
 *
 * @param node_forces Global nodal accumulator with at least XYZ components.
 * @param local_force Element force in node-major XYZ ordering.
 */
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

/**
 * Evaluates the reduced-integration continuum response and adds hourglass control.
 *
 * The common solid evaluation supplies all requested continuum quantities. The
 * reference hourglass operator is linear and state-neutral, so it contributes to
 * the complete tangent and to the internal force at the requested displacement,
 * but never to the separate perturbation geometric stiffness.
 */
MapMatrix C3D8R::evaluate(
    Precision*   tangent,
    Precision*   geometric_tangent,
    NodeData*    internal_force,
    const Field* displacement,
    const Field* linearization,
    bool         update_state
) {
    MapMatrix mapped = SolidElement<N>::evaluate(
        tangent,
        geometric_tangent,
        internal_force,
        displacement,
        linearization,
        update_state
    );

    const bool need_hourglass = tangent != nullptr || internal_force != nullptr;
    if (!need_hourglass) {
        return mapped;
    }

    const Matrix24 hourglass = hourglass_stiffness();

    if (internal_force != nullptr) {
        logging::error(displacement != nullptr,
            "C3D8R: hourglass force requires displacement");
        assemble_local_force(*internal_force, hourglass * local_displacement(*displacement));
    }

    if (tangent != nullptr) {
        MapMatrix full_tangent(tangent, ndof, ndof);
        full_tangent += hourglass;
        return full_tangent;
    }

    return mapped;
}

/**
 * @brief Returns the cached natural-space integration-point-to-node recovery map.
 *
 * Rows of the operator correspond to node_coords_local() order; columns follow
 * stress_strain_ip_rst() and therefore the constitutive quadrature/state order.
 * For one scalar component q, nodal recovery is q_node = E * q_ip. Vector and
 * tensor components use the same operator independently in their existing basis.
 *
 * A constant basis replicates the center-point value at all eight topology nodes.
 * math::extrapolate() fits coefficients through the source normal equations and
 * evaluates the polynomial at the natural node locations. Point-count and solver
 * validation are performed by that helper; no fallback is added here.
 *
 * Function-local static initialization builds the operator once from the fixed
 * natural topology and constitutive rule. No physical geometry, material history
 * or element runtime state is stored or changed by the recovery map.
 *
 * @return Shared immutable node-by-integration-point matrix in reference space.
 */
const RowMatrix& C3D8R::extrapolation_matrix() {
    // Fit the topology-selected polynomial basis at constitutive points and evaluate
    // it at the natural node coordinates. The static operator is reused by all instances.
    static const RowMatrix matrix = math::extrapolate(
        this->stress_strain_ip_rst(), this->node_coords_local(),
        {math::ExtrapolationBasis::F1});
    // Return the cached node-by-integration-point map without changing element state.
    return matrix;
}

} // namespace fem::model
