/**
 * @file c3d8r.cpp
 * @brief Implements the reduced-integration C3D8 solid and hourglass stabilization.
 *
 * The one-point continuum contribution uses the common solid material-point
 * state. Hourglass stiffness is obtained from a zero-strain constitutive tangent
 * without changing that physical history state, and the same hourglass matrix
 * supplies the residual and tangent contributions.
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#include "c3d8r.h"

#include <cmath>
#include <vector>

namespace fem::model {

C3D8R::C3D8R(ID elem_id, const std::array<ID, N>& node_ids)
    : C3D8(elem_id, node_ids) {}

std::string C3D8R::type_name() const {
    return "C3D8R";
}

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
 *     D_bar = (1 / V0) integral_A0 D dV0.
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

/**
 * Evaluates the zero-strain constitutive shear scale without advancing the
 * physical material-point history.
 *
 * The hourglass modulus is an auxiliary stabilization quantity rather than a
 * constitutive update of the current continuum state. The complete active state
 * row is therefore saved, the zero-strain tangent is evaluated through the same
 * material API, and the row is restored immediately afterwards.
 *
 * @return Mean material shear diagonal `(C44 + C55 + C66) / 3`.
 */
Precision C3D8R::hourglass_material_scale() {
    // Preserve the complete physical center-point history around the auxiliary
    // zero-strain tangent evaluation
    const Index state_row = this->mp_index(0);
    Field&      state_field = *this->_model_data->material_state;

    std::vector<Precision> saved_state(static_cast<std::size_t>(state_field.components));
    for (Index component = 0; component < state_field.components; ++component) {
        saved_state[static_cast<std::size_t>(component)] = state_field(state_row, component);
    }

    Precision* state = &state_field(state_row, 0);
    const Mat6 material_tangent = material_tangent_reference(
        Precision(0), Precision(0), Precision(0), state);

    for (Index component = 0; component < state_field.components; ++component) {
        state_field(state_row, component) = saved_state[static_cast<std::size_t>(component)];
    }

    const Precision shear_scale =
        (material_tangent(3, 3) + material_tangent(4, 4) + material_tangent(5, 5)) / Precision(3);

    logging::error(std::isfinite(shear_scale) && shear_scale > Precision(0),
        "C3D8R: invalid initial mean shear stiffness in element ", elem_id,
        "\nscale: ", shear_scale);

    return shear_scale;
}

/**
 * Builds the constant reference hourglass stabilization tangent.
 *
 * Primitive modes are projected with `I - D_bar X^T` so affine displacement
 * fields remain unstabilized. Their scalar stiffness uses the initial mean
 * constitutive shear diagonal, reference volume and mean-gradient norm:
 *
 *     k_hg = alpha G_eff V0 sum_a ||grad_bar N_a||^2.
 *
 * The scalar nodal matrix `k_hg G G^T` is expanded independently into all three
 * translational directions. The final matrix is symmetrized only to remove
 * round-off asymmetry.
 *
 * @return Constant 24-by-24 hourglass tangent in node-major DOF ordering.
 */
C3D8R::Matrix24 C3D8R::hourglass_stiffness() {
    const auto reference_coords = node_coords_reference();

    Precision reference_volume = Precision(0);
    const GradientMatrix mean_gradient = mean_reference_gradient(reference_volume);

    // Flanagan-Belytschko projection against the affine coordinate field:
    //
    //     G = (I - D_bar X^T) gamma.
    const StaticMatrix<N, N> projector =
        StaticMatrix<N, N>::Identity() - mean_gradient * reference_coords.transpose();
    const HourglassModes modes = projector * primitive_hourglass_modes();

    // Reference stabilization scale
    //
    //     k_hg = alpha G_eff V0 sum_a ||grad_bar N_a||^2.
    const Precision material_scale  = hourglass_material_scale();
    const Precision gradient_scale  = mean_gradient.array().square().sum();
    const Precision hourglass_scale =
        default_hourglass_coefficient * material_scale * reference_volume * gradient_scale;

    logging::error(std::isfinite(hourglass_scale) && hourglass_scale > Precision(0),
        "C3D8R: invalid hourglass stiffness in element ", elem_id,
        "\nscale: ", hourglass_scale);

    const StaticMatrix<N, N> scalar_stiffness =
        hourglass_scale * modes * modes.transpose();

    Matrix24 stiffness = Matrix24::Zero();

    // Expand the scalar matrix as K_hg = kron(H, I3).
    for (Index node_a = 0; node_a < N; ++node_a) {
        for (Index node_b = 0; node_b < N; ++node_b) {
            for (Dim dof = 0; dof < D; ++dof) {
                stiffness(D * node_a + dof, D * node_b + dof) =
                    scalar_stiffness(node_a, node_b);
            }
        }
    }

    return Precision(0.5) * (stiffness + stiffness.transpose());
}

/**
 * Collects current element translations relative to the reference geometry.
 *
 * @return Twenty-four-component displacement vector in node-major XYZ ordering.
 */
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
 * Assembles the one-point material stiffness with reference hourglass
 * stabilization.
 *
 * The inherited solid stiffness uses this element's virtual one-point
 * quadrature and advances the center-point state once. The auxiliary hourglass
 * tangent evaluates its material scale state-neutrally, then adds the constant
 * stabilization matrix.
 *
 * @param buffer Caller-provided dense 24-by-24 storage.
 * @return Mapped symmetric continuum-plus-hourglass stiffness.
 */
MapMatrix C3D8R::stiffness(Precision* buffer) {
    // The inherited continuum stiffness uses this element's virtual one-point
    // integration rule. Add the state-neutral hourglass tangent afterwards.
    MapMatrix mapped{buffer, ndof, ndof};

    C3D8::stiffness(buffer);
    mapped += hourglass_stiffness();
    mapped  = Precision(0.5) * (mapped + mapped.transpose());

    return mapped;
}

/**
 * Assembles the one-point continuum tangent and matching hourglass contribution.
 *
 * The common solid tangent performs the physical center-point constitutive
 * update exactly once. The auxiliary hourglass tangent is evaluated
 * state-neutrally and its matching linear force is added to the continuum
 * internal force.
 *
 * @param buffer Caller-provided dense 24-by-24 tangent storage.
 * @param ip_stress_state Global center-point PK2 stress field to update.
 * @param nodal_forces Global nodal internal-force field to increment.
 * @param displacement Trial displacement defining the current configuration.
 * @return Mapped consistent continuum-plus-hourglass tangent.
 */
MapMatrix C3D8R::stiffness_tangent(Precision*   buffer,
                                   Field&       ip_stress_state,
                                   NodeData&    nodal_forces,
                                   const Field& displacement) {
    const Matrix24 hourglass = hourglass_stiffness();

    MapMatrix mapped = SolidElement<N>::stiffness_tangent(
        buffer, ip_stress_state, nodal_forces, displacement);

    mapped += hourglass;
    mapped  = Precision(0.5) * (mapped + mapped.transpose());

    assemble_local_force(nodal_forces, hourglass * local_displacement());
    return mapped;
}

/**
 * Recovers nonlinear internal force from stored continuum stress and the
 * current hourglass displacement.
 *
 * The inherited Total-Lagrangian residual uses stored center-point PK2 stress.
 * The stabilization contribution is the exact derivative-compatible force
 * `f_hg = K_hg u_e` from the constant reference hourglass tangent.
 *
 * @param node_forces Global nodal internal-force field to increment.
 * @param ip_stress Stored center-point second Piola-Kirchhoff stress field.
 */
void C3D8R::compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) {
    // Continuum residual from the one-point Total-Lagrange formulation.
    C3D8::compute_internal_force_nonlinear(node_forces, ip_stress);

    // Matching constant reference hourglass force f_hg = K_hg u_e.
    const Matrix24 hourglass = hourglass_stiffness();
    assemble_local_force(node_forces, hourglass * local_displacement());
}

} // namespace fem::model
