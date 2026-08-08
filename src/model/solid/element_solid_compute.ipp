/**
 * @file element_solid_compute.ipp
 * @brief Implements solid stress recovery and nonlinear force/tangent assembly.
 *
 * Constitutive evaluations use the globally enumerated material-point state
 * associated with each solid integration point. The nonlinear tangent path
 * evaluates stress and material tangent together so an in-place history state
 * is advanced exactly once per material point and solver evaluation.
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "../../cos/rectangular_system.h"

namespace fem::model {

/**
 * Recovers solid strain and Cauchy stress at requested natural coordinates.
 *
 * Linearized recovery uses the reference-configuration small-strain B matrix.
 * Nonlinear recovery derives Green-Lagrange strain from the deformation
 * gradient, evaluates PK2 stress and pushes it forward to Cauchy stress for
 * output. Because requested coordinates may be nodal extrapolation points
 * rather than constitutive integration points, each location reuses the state
 * row of the nearest stiffness quadrature point in natural coordinates.
 *
 * Degenerate output mappings are skipped in the linearized path. At least one
 * of the optional output fields must be supplied.
 *
 * @param strain Optional output field for linearized or Green-Lagrange strain.
 * @param stress Optional output field for Cauchy stress.
 * @param displacement Global nodal displacement field.
 * @param rst Natural output coordinates, one point per row.
 * @param offset First output row belonging to this element.
 * @param use_green_lagrange_nl Select nonlinear Total-Lagrangian recovery.
 */
template<Index N>
void SolidElement<N>::compute_stress_strain(Field* strain,
                                            Field* stress,
                                            const Field& displacement,
                                            const RowMatrix& rst,
                                            int offset,
                                            bool use_green_lagrange_nl) {
    logging::error(strain != nullptr || stress != nullptr,
        "SolidElement: compute_stress_strain requires at least one output field");
    logging::error(rst.cols() >= 3,
        "SolidElement: stress/strain evaluation coordinates require at least 3 columns");

    auto reference_coords       = this->node_coords_reference();
    auto current_coords         = this->node_coords_current();
    auto local_displacement     = this->nodal_data<3>(displacement);
    auto local_disp_mat         = StaticMatrix<3, N>(local_displacement.transpose());
    auto local_displacement_vec = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);
    const auto& state_scheme    = this->integration_scheme_stiffness();

    for (Eigen::Index n = 0; n < rst.rows(); ++n) {
        const Precision r   = rst(n, 0);
        const Precision s   = rst(n, 1);
        const Precision t   = rst(n, 2);
        const Index     row = static_cast<Index>(offset + n);

        // Output coordinates may be nodal rather than constitutive integration
        // points. Reuse the history row of the nearest material point.
        Index state_ip = 0;
        auto  state_point = state_scheme.get_point(0);
        Precision state_distance =
            (r - state_point.r) * (r - state_point.r)
            + (s - state_point.s) * (s - state_point.s)
            + (t - state_point.t) * (t - state_point.t);

        for (Index ip = 1; ip < state_scheme.count(); ++ip) {
            state_point = state_scheme.get_point(ip);
            const Precision distance =
                (r - state_point.r) * (r - state_point.r)
                + (s - state_point.s) * (s - state_point.s)
                + (t - state_point.t) * (t - state_point.t);

            if (distance < state_distance) {
                state_ip       = ip;
                state_distance = distance;
            }
        }

        Precision* state = &(*this->_model_data->material_state)(this->mp_index(state_ip), 0);

        // Linearized recovery evaluates Cauchy stress directly from the
        // infinitesimal reference-configuration strain field
        if (!use_green_lagrange_nl) {
            Precision det;
            StaticMatrix<n_strain, D * N> B =
                this->strain_displacements(reference_coords, r, s, t, det, false);

            if (det <= Precision(0) || std::isnan(det) || std::isinf(det)) {
                continue;
            }

            const Vec6                   global_strain_voigt = B * local_displacement_vec;
            const VolumeStrainLinearized global_strain(global_strain_voigt);
            VolumeStressCauchy           global_stress;
            Mat6                         global_tangent;
            evaluate_material(r, s, t, global_strain, state, global_stress, global_tangent);

            for (Dim j = 0; j < n_strain; ++j) {
                if (strain) (*strain)(row, j) = global_strain.voigt()(j);
                if (stress) (*stress)(row, j) = global_stress.voigt()(j);
            }
            continue;
        }

        // Nonlinear recovery evaluates PK2 stress in the reference configuration
        // and pushes it forward for physical Cauchy-stress output
        const Mat3 F = this->deformation_gradient(reference_coords, current_coords, r, s, t);
        const VolumeStrainGreenLagrange green_lagrange =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);

        VolumeStressPK2 second_pk;
        Mat6            tangent;
        evaluate_material(r, s, t, green_lagrange, state, second_pk, tangent);

        const VolumeStressCauchy cauchy = second_pk.to_cauchy(F);

        for (Dim j = 0; j < n_strain; ++j) {
            if (strain) (*strain)(row, j) = green_lagrange.voigt()(j);
            if (stress) (*stress)(row, j) = cauchy.voigt()(j);
        }
    }
}

/**
 * Evaluates and stores the constitutive stress measure required by element
 * assembly at every solid integration point.
 *
 * Linearized evaluation stores Cauchy stress. The nonlinear Total-Lagrangian
 * path stores second Piola-Kirchhoff stress without a Cauchy push-forward because
 * subsequent internal-force and geometric-stiffness recovery require PK2. Each
 * integration point uses its own globally enumerated material-state row.
 *
 * @param stress_state Global integration-point stress field to update.
 * @param displacement Global nodal displacement field.
 * @param offset First output stress row supplied by the caller.
 * @param use_green_lagrange_nl Select PK2 or linearized Cauchy evaluation.
 */
template<Index N>
void SolidElement<N>::compute_stress_state(Field& stress_state,
                                           const Field& displacement,
                                           int offset,
                                           bool use_green_lagrange_nl) {
    const RowMatrix rst = stress_strain_ip_rst();
    if (rst.rows() == 0) {
        return;
    }

    // Prepare element geometry and the displacement vector once for all points
    const auto reference_coords       = this->node_coords_reference();
    const auto current_coords         = this->node_coords_current();
    const auto local_displacement     = this->nodal_data<3>(displacement);
    const auto local_disp_mat         = StaticMatrix<3, N>(local_displacement.transpose());
    const auto local_displacement_vec =
        Eigen::Map<const StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    for (Eigen::Index n = 0; n < rst.rows(); ++n) {
        const Precision r   = rst(n, 0);
        const Precision s   = rst(n, 1);
        const Precision t   = rst(n, 2);
        const Index     row = static_cast<Index>(offset + n);
        Precision*      state = &(*this->_model_data->material_state)(this->mp_index(static_cast<Index>(n)), 0);

        // Store linearized Cauchy stress directly when small-strain recovery is requested
        if (!use_green_lagrange_nl) {
            Precision det = Precision(0);
            const StaticMatrix<n_strain, D * N> B =
                this->strain_displacements(reference_coords, r, s, t, det, false);

            if (det <= Precision(0) || !std::isfinite(det)) {
                continue;
            }

            const Vec6                   strain_voigt = B * local_displacement_vec;
            const VolumeStrainLinearized strain(strain_voigt);
            VolumeStressCauchy           cauchy;
            Mat6                         tangent;
            evaluate_material(r, s, t, strain, state, cauchy, tangent);

            for (Dim component = 0; component < n_strain; ++component) {
                stress_state(row, component) = cauchy.voigt()(component);
            }
            continue;
        }

        // Store PK2 stress for the Total-Lagrangian residual and geometric tangent
        const Mat3 F = this->deformation_gradient(reference_coords, current_coords, r, s, t);
        const VolumeStrainGreenLagrange green_lagrange =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);

        VolumeStressPK2 second_pk;
        Mat6            tangent;
        evaluate_material(r, s, t, green_lagrange, state, second_pk, tangent);

        // Total-Lagrange internal forces and geometric stiffness require PK2.
        for (Dim component = 0; component < n_strain; ++component) {
            stress_state(row, component) = second_pk.voigt()(component);
        }
    }
}

/**
 * Assembles material tangent, geometric tangent and internal force from one
 * constitutive update at every solid material point.
 *
 * For each integration point the Green-Lagrange strain and PK2 stress satisfy
 * the Total-Lagrangian virtual-work relation. The material contribution is
 *
 *     K_mat = integral B^T C B dV0,
 *
 * while the stress-dependent geometric contribution is assembled from
 * `grad(N_a)^T S grad(N_b)`. The same PK2 stress is used for the internal force
 * `B^T S`, so no second constitutive evaluation of the material point is
 * required within the tangent assembly.
 *
 * @param buffer Caller-provided dense tangent storage.
 * @param ip_stress_state Global integration-point PK2 stress field to update.
 * @param nodal_forces Global nodal internal-force field to increment.
 * @param displacement Trial displacement defining the current configuration.
 * @return Complete material plus geometric tangent.
 */
template<Index N>
MapMatrix SolidElement<N>::stiffness_tangent(Precision*   buffer,
                                             Field&       ip_stress_state,
                                             NodeData&    nodal_forces,
                                             const Field& displacement) {
    (void) displacement;

    logging::error(ip_stress_state.components >= n_strain,
        "SolidElement: nonlinear stress state requires at least six components");
    logging::error(nodal_forces.components >= D,
        "SolidElement: nonlinear internal force requires at least three nodal components");

    const auto reference_coords = this->node_coords_reference();
    const auto current_coords   = this->node_coords_current();
    const auto& scheme          = this->integration_scheme_stiffness();

    StaticMatrix<D * N, D * N> tangent = StaticMatrix<D * N, D * N>::Zero();

    // Every quadrature point performs one constitutive update and immediately
    // uses the resulting stress and tangent for all element contributions.
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        Precision det0;
        const StaticMatrix<N, D> dN_dX = this->shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const Mat3 F = this->deformation_gradient(
            reference_coords, current_coords, point.r, point.s, point.t);
        const StaticMatrix<n_strain, D * N> B =
            this->green_lagrange_strain_displacement(dN_dX, F);
        const VolumeStrainGreenLagrange strain =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);

        Precision* state = &(*this->_model_data->material_state)(this->mp_index(ip), 0);

        VolumeStressPK2 stress;
        Mat6            material_tangent;
        evaluate_material(point.r, point.s, point.t, strain, state, stress, material_tangent);

        // Store the freshly evaluated PK2 stress in the common integration-point field.
        const Index stress_row = this->ip_index(ip);
        for (Index component = 0; component < n_strain; ++component) {
            ip_stress_state(stress_row, component) = stress.voigt()(component);
        }

        const Precision measure = point.w * det0;

        // Material tangent contribution B^T C B.
        tangent.noalias() += measure * B.transpose() * material_tangent * B;

        // Total-Lagrangian geometric tangent generated by the current PK2 stress.
        const Mat3 S = stress.tensor();
        for (Index a = 0; a < N; ++a) {
            const Vec3 dNa = dN_dX.row(a).transpose();

            for (Index b = 0; b < N; ++b) {
                const Vec3 dNb = dN_dX.row(b).transpose();
                const Precision s_ab = dNa.dot(S * dNb) * measure;

                for (Dim d = 0; d < D; ++d) {
                    tangent(D * a + d, D * b + d) += s_ab;
                }
            }
        }

        // Matching internal force contribution from the same constitutive state.
        const StaticVector<D * N> local_force = B.transpose() * stress.voigt() * measure;
        for (Index a = 0; a < N; ++a) {
            const Index node_id = static_cast<Index>(node_ids[a]);
            for (Dim d = 0; d < D; ++d) {
                nodal_forces(node_id, d) += local_force(D * a + d);
            }
        }
    }

    tangent = Precision(0.5) * (tangent + tangent.transpose());

    MapMatrix mapped{buffer, D * N, D * N};
    mapped = tangent;
    return mapped;
}

/**
 * Assembles the Total-Lagrangian internal force from a stored PK2 stress field.
 *
 * At each stiffness quadrature point the current deformation gradient builds
 * the nonlinear strain-displacement matrix. The virtual-work contribution
 * `B^T S dV0` is evaluated from the already stored stress so this recovery does
 * not call the material model or advance material history.
 *
 * @param node_forces Global nodal internal-force field to increment.
 * @param ip_stress Global integration-point PK2 stress field.
 */
template<Index N>
void SolidElement<N>::compute_internal_force_nonlinear(Field& node_forces,
                                                       const Field& ip_stress) {
    // Collect the two configurations and integration rule once for all points
    auto reference_coords = this->node_coords_reference();
    auto current_coords   = this->node_coords_current();
    auto scheme           = this->integration_scheme_stiffness();

    for (Index n = 0; n < scheme.count(); n++) {
        const Precision r = scheme.get_point(n).r;
        const Precision s = scheme.get_point(n).s;
        const Precision t = scheme.get_point(n).t;
        const Precision w = scheme.get_point(n).w;
        const Index row = this->ip_index(n);

        // Reconstruct the Total-Lagrangian B matrix in the reference configuration
        Precision det0;
        const auto dN_dX       = this->shape_derivatives_reference(reference_coords, r, s, t, det0);
        const auto F           = this->deformation_gradient(reference_coords, current_coords, r, s, t);
        const auto B           = this->green_lagrange_strain_displacement(dN_dX, F);

        // Integrate B^T S from the stored PK2 field without another material update
        const auto S_voigt     = ip_stress.row_vec6(row);
        const auto local_force = B.transpose() * S_voigt * (det0 * w);

        // Scatter the element contribution into global translational DOFs
        for (Index a = 0; a < N; ++a) {
            const Index node_id = static_cast<Index>(node_ids[a]);
            for (Index d = 0; d < D; ++d) {
                node_forces(node_id, d) += local_force(D * a + d);
            }
        }
    }
}

template<Index N>
void SolidElement<N>::compute_compliance(Field& displacement, Field& result) {
    Precision buffer[D * N * D * N];
    auto K = stiffness(buffer);

    auto local_disp_mat = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    Precision strain_energy = local_displacement.dot((K * local_displacement));
    result(elem_id, 0) = strain_energy;
}

/**
 * Computes the compliance derivative with respect to the three additional
 * material-orientation angles from the derivative of the constitutive tangent.
 *
 * With equilibrium `K u = f` and compliance `J = f^T u`, differentiation gives
 *
 *     J' = -u^T K' u
 *        = -integral eps^T C'_tan eps dV.
 *
 * @param displacement Current nodal displacement field.
 * @param result Element result field receiving the three angle derivatives.
 */
template<Index N>
void SolidElement<N>::compute_compliance_angle_derivative(Field& displacement, Field& result) {
    if (!this->_model_data || !this->_model_data->material_orientation) {
        return;
    }

    auto angles_field = this->_model_data->material_orientation;
    logging::error(angles_field->components == 3,
        "Field '", angles_field->name, "': material orientation requires 3 components");

    const Index row    = static_cast<Index>(this->elem_id);
    const Vec3  angles = angles_field->row_vec3(row);

    const Mat3 additional_rotation = cos::RectangularSystem::euler(
        angles(0), angles(1), angles(2)).get_axes(Vec3::Zero());

    const std::array<Mat3, 3> additional_rotation_derivatives {
        cos::RectangularSystem::derivative_rot_x(angles(0), angles(1), angles(2)),
        cos::RectangularSystem::derivative_rot_y(angles(0), angles(1), angles(2)),
        cos::RectangularSystem::derivative_rot_z(angles(0), angles(1), angles(2))
    };

    // Collect the element displacement and both configurations once
    auto local_disp_mat     = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    const auto reference_coords = this->node_coords_reference();
    const auto current_coords   = this->node_coords_current();
    const auto& scheme          = this->integration_scheme_stiffness();

    const Precision scaling    = element_stiffness_scale();
    Vec3            derivative = Vec3::Zero();

    // Integrate the energy sensitivity for all three rotation parameters
    for (Index n = 0; n < scheme.count(); ++n) {
        const auto point = scheme.get_point(n);
        Precision det;

        const StaticMatrix<n_strain, D * N> B =
            this->strain_displacements(current_coords, point.r, point.s, point.t, det);
        const StaticVector<n_strain> strain = B * local_displacement;

        // Evaluate transformation derivatives at the physical reference point
        // using the material state row of this integration point
        const Vec3 position_reference =
            this->interpolate<D>(reference_coords, point.r, point.s, point.t);
        Precision* state = &(*this->_model_data->material_state)(this->mp_index(n), 0);

        const auto tangent_derivatives = get_section()->tangent_rotation_derivatives(
            position_reference,
            additional_rotation,
            additional_rotation_derivatives,
            state
        );

        // Accumulate epsilon^T C'_tan epsilon with the physical volume measure
        for (Index i = 0; i < 3; ++i) {
            derivative(i) += scaling * point.w * strain.dot(tangent_derivatives[i] * strain) * det;
        }
    }

    result(elem_id, 0) = derivative(0);
    result(elem_id, 1) = derivative(1);
    result(elem_id, 2) = derivative(2);
}

} // namespace fem::model
