/**
 * @file c3d8r.cpp
 * @brief Implements the reduced-integration C3D8 solid and physical hourglass stabilization.
 *
 * The one-point continuum contribution uses the common solid material-point
 * state. Hourglass stabilization restores only the missing deviatoric reference
 * stiffness inside the twelve-dimensional hourglass subspace.
 *
 * @author Finn Eggers
 * @date 26.09.2026
 */

#include "c3d8r.h"

#include <Eigen/LU>

#include <cmath>

namespace fem::model {

C3D8R::C3D8R(ID elem_id, const std::array<ID, N>& node_ids)
    : C3D8(elem_id, node_ids) {}

std::string C3D8R::type_name() const {
    return "C3D8R";
}

const math::quadrature::Quadrature& C3D8R::integration_scheme_stiffness() const {
    static const math::quadrature::Quadrature quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_CONSTANT
    };
    return quadrature;
}

RowMatrix C3D8R::stress_strain_nodal_rst() {
    return RowMatrix::Zero(N, D);
}

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
 * The average gradient is used only to project the primitive hourglass patterns
 * away from affine coordinate fields. The physical stabilization metric itself
 * is assembled separately from full and reduced reference integration.
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
            reference_coords, point.r, point.s, point.t, det0);

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
 * Builds the twelve-dimensional translational hourglass basis.
 *
 * The four scalar Flanagan-Belytschko modes are projected with
 *
 *     G = (I - D_bar X^T) gamma
 *
 * and embedded independently in x, y and z displacement directions.
 */
C3D8R::HourglassBasis C3D8R::hourglass_basis() {
    const auto reference_coords = node_coords_reference();
    const GradientMatrix mean_gradient = mean_reference_gradient();

    const StaticMatrix<N, N> projector =
        StaticMatrix<N, N>::Identity() - mean_gradient * reference_coords.transpose();
    const HourglassModes modes = projector * primitive_hourglass_modes();

    HourglassBasis basis = HourglassBasis::Zero();

    for (Index mode = 0; mode < 4; ++mode) {
        for (Dim dof = 0; dof < D; ++dof) {
            const Index column = D * mode + dof;
            for (Index node = 0; node < N; ++node) {
                basis(D * node + dof, column) = modes(node, mode);
            }
        }
    }

    return basis;
}

/**
 * Returns the generic constitutive tangent restricted to deviatoric strain and
 * stress subspaces.
 *
 * No material parameters are inspected by C3D8R. The ordinary material API
 * supplies the zero-strain reference tangent, after which a kinematic
 * deviatoric projector removes the hydrostatic strain/stress direction:
 *
 *     C_hg = P_dev^T C P_dev.
 *
 * This keeps nearly incompressible bulk stiffness out of the stabilization
 * without special-casing isotropic, hyperelastic or plastic material models.
 */
Mat6 C3D8R::deviatoric_reference_tangent() {
    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);

    const Mat6 material_tangent = material_tangent_reference(
        Precision(0), Precision(0), Precision(0), old_state, nullptr);

    Mat6 projector = Mat6::Identity();
    for (Index row = 0; row < 3; ++row) {
        for (Index column = 0; column < 3; ++column) {
            projector(row, column) -= Precision(1) / Precision(3);
        }
    }

    Mat6 tangent = projector.transpose() * material_tangent * projector;
    tangent = Precision(0.5) * (tangent + tangent.transpose());

    logging::error(tangent.allFinite(),
        "C3D8R: invalid deviatoric reference tangent in element ", elem_id);

    return tangent;
}

/**
 * Builds the reference hourglass stabilization tangent.
 *
 * A full 2x2x2 rule and the actual one-point rule are evaluated with the same
 * generic deviatoric constitutive tangent. Their difference is the stiffness
 * omitted by reduced integration:
 *
 *     K_missing = K_dev,8GP - K_dev,1GP.
 *
 * Only the component of this operator acting inside the twelve-dimensional
 * hourglass subspace is retained. For a non-orthonormal basis H,
 *
 *     M    = H^T H,
 *     K_q  = H^T K_missing H,
 *     K_hg = H M^-1 K_q M^-1 H^T.
 *
 * This replaces the previous scalar alpha*G*geometry metric and introduces no
 * empirical hourglass coefficient.
 */
C3D8R::Matrix24 C3D8R::hourglass_stiffness() {
    const auto reference_coords = node_coords_reference();
    const HourglassBasis basis = hourglass_basis();
    const Mat6 material_tangent = deviatoric_reference_tangent();

    Matrix24 full_stiffness    = Matrix24::Zero();
    Matrix24 reduced_stiffness = Matrix24::Zero();

    static const math::quadrature::Quadrature full_quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };

    for (Index q = 0; q < full_quadrature.count(); ++q) {
        const auto point = full_quadrature.get_point(q);

        Precision det0 = Precision(0);
        const auto gradient = shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const StaticMatrix<6, ndof> B = strain_displacement(gradient);

        logging::error(std::isfinite(det0) && det0 > Precision(0),
            "C3D8R: invalid full-integration reference determinant in element ", elem_id,
            "\ndet(J0): ", det0);

        full_stiffness.noalias() +=
            B.transpose() * material_tangent * B * (det0 * point.w);
    }

    const auto& reduced_quadrature = integration_scheme_stiffness();
    for (Index q = 0; q < reduced_quadrature.count(); ++q) {
        const auto point = reduced_quadrature.get_point(q);

        Precision det0 = Precision(0);
        const auto gradient = shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const StaticMatrix<6, ndof> B = strain_displacement(gradient);

        logging::error(std::isfinite(det0) && det0 > Precision(0),
            "C3D8R: invalid reduced-integration reference determinant in element ", elem_id,
            "\ndet(J0): ", det0);

        reduced_stiffness.noalias() +=
            B.transpose() * material_tangent * B * (det0 * point.w);
    }

    const Matrix24 missing_stiffness =
        Precision(0.5) * ((full_stiffness - reduced_stiffness)
                        + (full_stiffness - reduced_stiffness).transpose());

    const Matrix12 gram = basis.transpose() * basis;
    const Matrix12 gram_inverse = gram.inverse();

    logging::error(gram_inverse.allFinite(),
        "C3D8R: singular hourglass basis in element ", elem_id);

    Matrix12 modal_stiffness = basis.transpose() * missing_stiffness * basis;
    modal_stiffness = Precision(0.5) * (modal_stiffness + modal_stiffness.transpose());

    Matrix24 stiffness =
        basis * gram_inverse * modal_stiffness * gram_inverse * basis.transpose();
    stiffness = Precision(0.5) * (stiffness + stiffness.transpose());

    logging::error(stiffness.allFinite(),
        "C3D8R: invalid hourglass stiffness in element ", elem_id);

    return stiffness;
}

C3D8R::Vector24 C3D8R::local_displacement(const Field& displacement) {
    const GradientMatrix local = this->nodal_data<D>(displacement);
    Vector24 result = Vector24::Zero();

    for (Index node = 0; node < N; ++node) {
        for (Dim dof = 0; dof < D; ++dof) {
            result(D * node + dof) = local(node, dof);
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

MapMatrix C3D8R::stiffness(Precision* buffer) {
    MapMatrix mapped{buffer, ndof, ndof};

    C3D8::stiffness(buffer);
    mapped += hourglass_stiffness();
    mapped  = Precision(0.5) * (mapped + mapped.transpose());

    return mapped;
}

MapMatrix C3D8R::stiffness_tangent(Precision*   buffer,
                                   NodeData&    nodal_forces,
                                   const Field& displacement) {
    const Matrix24 hourglass = hourglass_stiffness();

    MapMatrix mapped = SolidElement<N>::stiffness_tangent(buffer, nodal_forces, displacement);

    assemble_local_force(nodal_forces, hourglass * local_displacement(displacement));

    if (buffer == nullptr) {
        return mapped;
    }

    mapped += hourglass;
    return mapped;
}

} // namespace fem::model
