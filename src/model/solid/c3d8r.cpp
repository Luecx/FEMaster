/**
 * @file c3d8r.cpp
 * @brief Implements the reduced-integration C3D8 solid and projected hourglass stabilization.
 *
 * The one-point continuum contribution uses the common solid material-point
 * state. Hourglass stiffness is obtained from the fully integrated deviatoric
 * zero-strain reference tangent projected onto the twelve non-affine element
 * modes. The resulting constant matrix supplies matching residual and tangent
 * contributions without an empirical stabilization coefficient.
 *
 * @author Finn Eggers
 * @date 13.08.2026
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
 * The mean gradient is used only to remove physical affine displacement fields
 * from the primitive hourglass basis.
 *
 * @return Mean reference gradient with one row per element node.
 */
C3D8R::GradientMatrix C3D8R::mean_reference_gradient() {
    const auto reference_coords = node_coords_reference();

    static const math::quadrature::Quadrature full_quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };

    GradientMatrix integrated_gradient = GradientMatrix::Zero();
    Precision reference_volume = Precision(0);

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
 * Builds the parameter-free reference hourglass stabilization tangent.
 *
 * The projected Flanagan-Belytschko modes span the four scalar non-affine nodal
 * patterns. Their orthogonal scalar projector is expanded independently into the
 * three displacement directions, giving the twelve-dimensional projector
 * `P_hg`.
 *
 * A full two-by-two-by-two rule then integrates the zero-strain reference
 * stiffness using only the deviatoric part of the constitutive tangent:
 *
 *     K_dev^(8GP) = integral B^T P_dev^T C0 P_dev B dV0.
 *
 * Finally,
 *
 *     K_hg = P_hg K_dev^(8GP) P_hg.
 *
 * Consequently no bulk stiffness and no empirical hourglass coefficient enter
 * the stabilization. `K_hg` is constant in the reference configuration and is
 * cached after the first evaluation.
 *
 * @return Constant 24-by-24 projected hourglass tangent.
 */
C3D8R::Matrix24 C3D8R::hourglass_stiffness() {
    if (hourglass_stiffness_cached) {
        return hourglass_stiffness_cache;
    }

    const auto reference_coords = node_coords_reference();
    const GradientMatrix mean_gradient = mean_reference_gradient();

    // Remove physical affine scalar fields from the primitive modes.
    const StaticMatrix<N, N> affine_projector =
        StaticMatrix<N, N>::Identity() - mean_gradient * reference_coords.transpose();
    const HourglassModes modes = affine_projector * primitive_hourglass_modes();

    // Orthogonal projector onto the four scalar non-affine modes:
    //
    //     P_s = G (G^T G)^-1 G^T.
    const StaticMatrix<4, 4> gram = modes.transpose() * modes;
    logging::error(std::isfinite(gram.determinant()) && std::abs(gram.determinant()) > Precision(1e-12),
        "C3D8R: singular projected hourglass basis in element ", elem_id,
        "\ndet(G^T G): ", gram.determinant());

    const StaticMatrix<N, N> scalar_projector =
        modes * gram.inverse() * modes.transpose();

    Matrix24 hourglass_projector = Matrix24::Zero();
    for (Index node_a = 0; node_a < N; ++node_a) {
        for (Index node_b = 0; node_b < N; ++node_b) {
            for (Dim dof = 0; dof < D; ++dof) {
                hourglass_projector(D * node_a + dof, D * node_b + dof) =
                    scalar_projector(node_a, node_b);
            }
        }
    }

    // Deviatoric projector in the solid engineering-Voigt convention
    // [11,22,33,23,13,12]. Normal strains are made trace-free while engineering
    // shear strains remain unchanged.
    Mat6 deviatoric_projector = Mat6::Zero();
    for (Dim i = 0; i < 3; ++i) {
        for (Dim j = 0; j < 3; ++j) {
            deviatoric_projector(i, j) =
                (i == j ? Precision(1) : Precision(0)) - Precision(1) / Precision(3);
        }
    }
    deviatoric_projector(3, 3) = Precision(1);
    deviatoric_projector(4, 4) = Precision(1);
    deviatoric_projector(5, 5) = Precision(1);

    // Preserve the physical center-point material state around all auxiliary
    // zero-strain tangent evaluations. Each quadrature point starts from exactly
    // the same saved state, and the original row is restored afterwards.
    const Index state_row = this->mp_index(0);
    Field& state_field = *this->_model_data->material_state;

    std::vector<Precision> saved_state(static_cast<std::size_t>(state_field.components));
    for (Index component = 0; component < state_field.components; ++component) {
        saved_state[static_cast<std::size_t>(component)] = state_field(state_row, component);
    }

    auto restore_state = [&]() {
        for (Index component = 0; component < state_field.components; ++component) {
            state_field(state_row, component) = saved_state[static_cast<std::size_t>(component)];
        }
    };

    Precision* state = &state_field(state_row, 0);

    static const math::quadrature::Quadrature full_quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };

    Matrix24 deviatoric_reference_stiffness = Matrix24::Zero();

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
        const auto B = strain_displacement(gradient);

        restore_state();
        const Mat6 material_tangent = material_tangent_reference(
            point.r,
            point.s,
            point.t,
            state
        );
        const Mat6 deviatoric_tangent =
            deviatoric_projector.transpose() * material_tangent * deviatoric_projector;

        const Precision measure = det0 * point.w;
        deviatoric_reference_stiffness +=
            B.transpose() * (deviatoric_tangent * B) * measure;
    }

    restore_state();

    hourglass_stiffness_cache =
        hourglass_projector * deviatoric_reference_stiffness * hourglass_projector;
    hourglass_stiffness_cache =
        Precision(0.5) * (hourglass_stiffness_cache + hourglass_stiffness_cache.transpose());
    hourglass_stiffness_cached = true;

    return hourglass_stiffness_cache;
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
 * Assembles the one-point material stiffness with projected reference hourglass
 * stabilization.
 *
 * @param buffer Caller-provided dense 24-by-24 storage.
 * @return Mapped symmetric continuum-plus-hourglass stiffness.
 */
MapMatrix C3D8R::stiffness(Precision* buffer) {
    MapMatrix mapped{buffer, ndof, ndof};

    C3D8::stiffness(buffer);
    mapped += hourglass_stiffness();
    mapped  = Precision(0.5) * (mapped + mapped.transpose());

    return mapped;
}

/**
 * Assembles the one-point continuum tangent and matching projected hourglass
 * contribution.
 *
 * The common solid tangent performs the physical center-point constitutive
 * update exactly once. The constant reference hourglass tangent adds its exact
 * matching linear force.
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
 * current projected hourglass displacement.
 *
 * @param node_forces Global nodal internal-force field to increment.
 * @param ip_stress Stored center-point second Piola-Kirchhoff stress field.
 */
void C3D8R::compute_internal_force_nonlinear(Field& node_forces, const Field& ip_stress) {
    C3D8::compute_internal_force_nonlinear(node_forces, ip_stress);

    const Matrix24 hourglass = hourglass_stiffness();
    assemble_local_force(node_forces, hourglass * local_displacement());
}

} // namespace fem::model
