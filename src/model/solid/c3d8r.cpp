/**
 * @file c3d8r.cpp
 * @brief Implements reduced C3D8 integration with objective finite-strain hourglass control.
 *
 * The physical material response remains a one-point C3D8R evaluation. The
 * stabilization samples only kinematics at the ordinary 2x2x2 C3D8 points and
 * penalizes both non-constant deviatoric Green-Lagrange strain and non-uniform
 * volume change relative to the element center. No additional material states
 * are created.
 *
 * @author Finn Eggers
 * @date 26.09.2026
 */

#include "c3d8r.h"

#include <cmath>

namespace fem::model {

namespace {

/**
 * Extracts one basis-invariant deviatoric stiffness scale from the generic
 * engineering-Voigt tangent. The five vectors form an orthonormal basis of the
 * symmetric deviatoric tensor space. Each isotropic basis energy equals 2G, so
 *
 *     mu_eff = (1/10) sum_i e_i^T C_dev e_i
 *
 * recovers exactly the shear modulus G without inspecting the material model.
 */
Precision effective_deviatoric_modulus(const Mat6& tangent) {
    const Precision inv_sqrt_2 = Precision(1) / std::sqrt(Precision(2));
    const Precision inv_sqrt_6 = Precision(1) / std::sqrt(Precision(6));
    const Precision sqrt_2     = std::sqrt(Precision(2));

    std::array<Vec6, 5> basis;
    basis[0] <<  inv_sqrt_2, -inv_sqrt_2, Precision(0), Precision(0), Precision(0), Precision(0);
    basis[1] <<  inv_sqrt_6,  inv_sqrt_6, -Precision(2) * inv_sqrt_6, Precision(0), Precision(0), Precision(0);
    basis[2] <<  Precision(0), Precision(0), Precision(0), sqrt_2, Precision(0), Precision(0);
    basis[3] <<  Precision(0), Precision(0), Precision(0), Precision(0), sqrt_2, Precision(0);
    basis[4] <<  Precision(0), Precision(0), Precision(0), Precision(0), Precision(0), sqrt_2;

    Precision sum = Precision(0);
    for (const auto& direction : basis) {
        sum += direction.dot(tangent * direction);
    }

    return sum / Precision(10);
}

} // namespace

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

/**
 * Returns a state-neutral constitutive tangent restricted to the deviatoric
 * strain/stress subspace.
 *
 * C3D8R never inspects concrete material parameters. The generic material
 * tangent is queried at zero Green-Lagrange strain from committed state and
 * projected kinematically:
 *
 *     C_hg = P_dev^T C P_dev.
 *
 * The hydrostatic direction is excluded so nearly incompressible bulk stiffness
 * cannot enter the hourglass control.
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
 * Builds the infinitesimal tangent of the finite hourglass formulation.
 *
 * At each full-integration point the non-constant strain operator is
 *
 *     A_q = B_q - B_0,
 *
 * so every affine displacement field remains exactly unstabilized. The linear
 * hourglass operator combines the deviatoric missing-strain term with the
 * linearization of the non-uniform logarithmic volume ratio:
 *
 *     K_hg = sum_q [A_q^T C_dev A_q
 *                  + mu_eff a_q a_q^T] dV0,
 *
 * where a_q = A_q^T [1,1,1,0,0,0]^T. This is the zero-deformation tangent of
 * finite_hourglass().
 */
C3D8R::Matrix24 C3D8R::linear_hourglass_stiffness() {
    const auto reference_coords = node_coords_reference();
    const Mat6 material_tangent = deviatoric_reference_tangent();
    const Precision shear_scale = effective_deviatoric_modulus(material_tangent);

    logging::error(std::isfinite(shear_scale) && shear_scale > Precision(0),
        "C3D8R: invalid deviatoric stiffness scale in element ", elem_id,
        "\nscale: ", shear_scale);

    Vec6 volume_direction = Vec6::Zero();
    volume_direction.head<3>().setOnes();

    Precision center_det = Precision(0);
    const auto center_gradient = shape_derivatives_reference(
        reference_coords, Precision(0), Precision(0), Precision(0), center_det);
    const StaticMatrix<6, ndof> center_B = strain_displacement(center_gradient);

    logging::error(std::isfinite(center_det) && center_det > Precision(0),
        "C3D8R: invalid center reference determinant in element ", elem_id,
        "\ndet(J0): ", center_det);

    static const math::quadrature::Quadrature full_quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };

    Matrix24 stiffness = Matrix24::Zero();

    for (Index q = 0; q < full_quadrature.count(); ++q) {
        const auto point = full_quadrature.get_point(q);

        Precision det0 = Precision(0);
        const auto gradient = shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const StaticMatrix<6, ndof> B = strain_displacement(gradient);
        const StaticMatrix<6, ndof> A = B - center_B;

        logging::error(std::isfinite(det0) && det0 > Precision(0),
            "C3D8R: invalid reference determinant in element ", elem_id,
            "\ndet(J0): ", det0);

        const Precision measure = det0 * point.w;
        const Vector24 volume_mode = A.transpose() * volume_direction;

        stiffness.noalias() +=
            (A.transpose() * material_tangent * A
             + shear_scale * volume_mode * volume_mode.transpose()) * measure;
    }

    stiffness = Precision(0.5) * (stiffness + stiffness.transpose());

    logging::error(stiffness.allFinite(),
        "C3D8R: invalid linear hourglass stiffness in element ", elem_id);

    return stiffness;
}

/**
 * Evaluates objective finite-strain hourglass force and optional tangent.
 *
 * The full 2x2x2 points are kinematic sampling points only; the physical
 * constitutive state remains the single reduced point at the element center.
 * For each sample,
 *
 *     E_hg,q = E_q - E_0,
 *     theta_q = log(J_q / J_0),
 *
 * and the stabilization energy is
 *
 *     Psi_hg = 1/2 sum_q [E_hg,q^T C_dev E_hg,q
 *                         + mu_eff theta_q^2] dV0.
 *
 * The first term controls non-constant deviatoric strain. The logarithmic
 * volume-ratio term is exactly zero for every homogeneous deformation but tends
 * to infinity as a local sample collapses while the center remains regular.
 * Its scale mu_eff is derived from C_dev and equals G for isotropic materials.
 *
 * Green-Lagrange strain and J make the stabilization objective under arbitrary
 * rigid rotation.
 *
 * When requested, the tangent contains both the material term and the exact
 * Total-Lagrangian geometric terms generated by the q and center strain
 * variations. The routine also evaluates F at all eight sample points; any
 * local inversion therefore throws before a folded C3D8R trial can be accepted.
 */
void C3D8R::finite_hourglass(const Field& displacement,
                             Vector24&    local_force,
                             Matrix24*    tangent) {
    const auto reference_coords   = node_coords_reference();
    const auto local_displacement = this->nodal_data<D>(displacement);
    const auto current_coords     = reference_coords + local_displacement;
    const Mat6 material_tangent   = deviatoric_reference_tangent();
    const Precision shear_scale   = effective_deviatoric_modulus(material_tangent);

    logging::error(std::isfinite(shear_scale) && shear_scale > Precision(0),
        "C3D8R: invalid deviatoric stiffness scale in element ", elem_id,
        "\nscale: ", shear_scale);

    Precision center_det0 = Precision(0);
    const auto center_gradient = shape_derivatives_reference(
        reference_coords, Precision(0), Precision(0), Precision(0), center_det0);
    const Mat3 center_F = deformation_gradient(
        reference_coords, current_coords, Precision(0), Precision(0), Precision(0));
    const Precision center_J = center_F.determinant();
    const Mat3 center_F_inv_T = center_F.inverse().transpose();
    const Vec6 center_strain =
        VolumeStrainGreenLagrange::from_deformation_gradient(center_F).voigt();
    const StaticMatrix<6, ndof> center_B =
        green_lagrange_strain_displacement(center_gradient, center_F);

    logging::error(std::isfinite(center_det0) && center_det0 > Precision(0),
        "C3D8R: invalid center reference determinant in element ", elem_id,
        "\ndet(J0): ", center_det0);

    static const math::quadrature::Quadrature full_quadrature{
        math::quadrature::DOMAIN_ISO_HEX,
        math::quadrature::ORDER_QUADRATIC
    };

    local_force.setZero();
    if (tangent != nullptr) {
        tangent->setZero();
    }

    for (Index q = 0; q < full_quadrature.count(); ++q) {
        const auto point = full_quadrature.get_point(q);

        Precision det0 = Precision(0);
        const auto gradient = shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const Mat3 F = deformation_gradient(
            reference_coords, current_coords, point.r, point.s, point.t);
        const Precision J = F.determinant();

        logging::error(std::isfinite(J) && J > Precision(0),
            "C3D8R: inverted full-integration sample in element ", elem_id,
            "\npoint: (", point.r, ", ", point.s, ", ", point.t, ")",
            "\ndet(F): ", J);

        const Vec6 strain =
            VolumeStrainGreenLagrange::from_deformation_gradient(F).voigt();
        const StaticMatrix<6, ndof> B =
            green_lagrange_strain_displacement(gradient, F);
        const StaticMatrix<6, ndof> A = B - center_B;

        const Vec6 hourglass_strain = strain - center_strain;
        const Vec6 hourglass_stress = material_tangent * hourglass_strain;
        const Precision theta = std::log(J / center_J);
        const Precision measure = det0 * point.w;

        const Mat3 F_inv_T = F.inverse().transpose();
        Vector24 volume_gradient = Vector24::Zero();

        std::array<Vec3, N> current_volume_gradients;
        std::array<Vec3, N> center_volume_gradients;

        for (Index a = 0; a < N; ++a) {
            const Vec3 dNa_q = gradient.row(a).transpose();
            const Vec3 dNa_0 = center_gradient.row(a).transpose();

            current_volume_gradients[a] = F_inv_T * dNa_q;
            center_volume_gradients[a]  = center_F_inv_T * dNa_0;

            for (Dim d = 0; d < D; ++d) {
                volume_gradient(D * a + d) =
                    current_volume_gradients[a](d) - center_volume_gradients[a](d);
            }
        }

        local_force.noalias() +=
            (A.transpose() * hourglass_stress
             + shear_scale * theta * volume_gradient) * measure;

        if (tangent == nullptr) {
            continue;
        }

        tangent->noalias() +=
            (A.transpose() * material_tangent * A
             + shear_scale * volume_gradient * volume_gradient.transpose()) * measure;

        const Mat3 S = VolumeStressPK2(hourglass_stress).tensor();

        // Exact second variation of E_q - E_0. The q contribution has positive
        // sign and the center contribution negative sign.
        for (Index a = 0; a < N; ++a) {
            const Vec3 dNa_q = gradient.row(a).transpose();
            const Vec3 dNa_0 = center_gradient.row(a).transpose();

            for (Index b = 0; b < N; ++b) {
                const Vec3 dNb_q = gradient.row(b).transpose();
                const Vec3 dNb_0 = center_gradient.row(b).transpose();
                const Precision coefficient =
                    (dNa_q.dot(S * dNb_q) - dNa_0.dot(S * dNb_0)) * measure;

                for (Dim d = 0; d < D; ++d) {
                    (*tangent)(D * a + d, D * b + d) += coefficient;
                }

                const Mat3 volume_hessian =
                    -current_volume_gradients[b] * current_volume_gradients[a].transpose()
                    +center_volume_gradients[b] * center_volume_gradients[a].transpose();

                tangent->block<3, 3>(D * a, D * b).noalias() +=
                    shear_scale * theta * volume_hessian * measure;
            }
        }
    }

    logging::error(local_force.allFinite(),
        "C3D8R: invalid finite hourglass force in element ", elem_id);

    if (tangent != nullptr) {
        *tangent = Precision(0.5) * (*tangent + tangent->transpose());
        logging::error(tangent->allFinite(),
            "C3D8R: invalid finite hourglass tangent in element ", elem_id);
    }
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
    mapped += linear_hourglass_stiffness();
    mapped  = Precision(0.5) * (mapped + mapped.transpose());

    return mapped;
}

MapMatrix C3D8R::stiffness_tangent(Precision*   buffer,
                                   NodeData&    nodal_forces,
                                   const Field& displacement) {
    Vector24 hourglass_force;
    Matrix24 hourglass_tangent;

    // Evaluate the full-point geometry and hourglass response first. A locally
    // inverted trial throws here and is handled by nonlinear load-control
    // cutback before the physical center material state is advanced.
    finite_hourglass(
        displacement,
        hourglass_force,
        buffer != nullptr ? &hourglass_tangent : nullptr
    );

    MapMatrix mapped =
        SolidElement<N>::stiffness_tangent(buffer, nodal_forces, displacement);

    assemble_local_force(nodal_forces, hourglass_force);

    if (buffer == nullptr) {
        return mapped;
    }

    mapped += hourglass_tangent;
    return mapped;
}

} // namespace fem::model
