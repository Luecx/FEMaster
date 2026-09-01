/**
 * @file element_solid_compute.ipp
 * @brief Implements solid stress recovery and nonlinear force/tangent assembly.
 *
 * Constitutive evaluations use globally enumerated committed/trial material-state
 * rows associated with each solid integration point. Result recovery and other
 * auxiliary evaluations are state-neutral. The physical nonlinear tangent path
 * evaluates PK2 stress and material tangent together exactly once per material
 * point and immediately consumes both quantities locally.
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "../../cos/rectangular_system.h"
#include "../../section/section_solid.h"

namespace fem::model {

/**
 * Recovers solid strain and Cauchy stress at requested natural coordinates.
 *
 * Linearized recovery uses the reference-configuration small-strain B matrix.
 * Nonlinear recovery derives Green-Lagrange strain from the deformation
 * gradient, evaluates PK2 stress and pushes it forward to Cauchy stress for
 * output. Because requested coordinates may be nodal extrapolation points
 * rather than constitutive integration points, each location reuses the
 * committed state row of the nearest stiffness quadrature point in natural
 * coordinates.
 *
 * Recovery is strictly state-neutral: constitutive history may be read but no
 * target trial state is supplied. Degenerate output mappings are skipped in the
 * linearized path. At least one optional output field must be supplied.
 *
 * @param strain Optional output field for linearized or Green-Lagrange strain.
 * @param stress Optional output field for Cauchy stress.
 * @param displacement Global nodal displacement field.
 * @param rst Natural output coordinates, one point per row.
 * @param offset First output row belonging to this element.
 * @param use_green_lagrange_nl Select nonlinear Total-Lagrangian recovery.
 */
template<Index N>
void SolidElement<N>::compute_stress_strain(Field*           strain,
                                            Field*           stress,
                                            const Field&     displacement,
                                            const RowMatrix& rst,
                                            int              offset,
                                            bool             use_green_lagrange_nl) {
    logging::error(strain != nullptr || stress != nullptr,
        "SolidElement: compute_stress_strain requires at least one output field");
    logging::error(rst.cols() >= 3,
        "SolidElement: stress/strain evaluation coordinates require at least 3 columns");

    const auto reference_coords       = this->node_coords_reference();
    const auto local_displacement     = this->nodal_data<3>(displacement);
    const auto local_disp_mat         = StaticMatrix<3, N>(local_displacement.transpose());
    const auto local_displacement_vec = Eigen::Map<const StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);
    const auto current_coords         = reference_coords + local_displacement;
    const auto& state_scheme          = this->integration_scheme_stiffness();

    for (Eigen::Index n = 0; n < rst.rows(); ++n) {
        const Precision r   = rst(n, 0);
        const Precision s   = rst(n, 1);
        const Precision t   = rst(n, 2);
        const Index     row = static_cast<Index>(offset + n);

        // Output coordinates may be nodal rather than constitutive integration
        // points. Reuse the committed history row of the nearest material point.
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

        const Index      state_row = this->mp_index(state_ip);
        const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);

        // Linearized recovery evaluates Cauchy stress directly from the
        // infinitesimal reference-configuration strain field.
        if (!use_green_lagrange_nl) {
            Precision det;
            const StaticMatrix<n_strain, D * N> B =
                this->strain_displacements(reference_coords, r, s, t, det, false);

            if (det <= Precision(0) || !std::isfinite(det)) {
                continue;
            }

            const Vec6                   global_strain_voigt = B * local_displacement_vec;
            const VolumeStrainLinearized global_strain(global_strain_voigt);
            VolumeStressCauchy           global_stress;
            Mat6                         global_tangent;
            evaluate_material(r, s, t, global_strain, old_state, nullptr, global_stress, global_tangent);

            for (Dim j = 0; j < n_strain; ++j) {
                if (strain) (*strain)(row, j) = global_strain.voigt()(j);
                if (stress) (*stress)(row, j) = global_stress.voigt()(j);
            }
            continue;
        }

        // Nonlinear recovery evaluates PK2 stress in the reference configuration
        // and pushes it forward for physical Cauchy-stress output.
        const Mat3 F = this->deformation_gradient(reference_coords, current_coords, r, s, t);
        const VolumeStrainGreenLagrange green_lagrange =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);

        VolumeStressPK2 second_pk;
        Mat6            tangent;
        evaluate_material(r, s, t, green_lagrange, old_state, nullptr, second_pk, tangent);

        const VolumeStressCauchy cauchy = second_pk.to_cauchy(F);

        for (Dim j = 0; j < n_strain; ++j) {
            if (strain) (*strain)(row, j) = green_lagrange.voigt()(j);
            if (stress) (*stress)(row, j) = cauchy.voigt()(j);
        }
    }
}

/**
 * Assembles Total-Lagrangian internal force and, when requested, the complete
 * consistent solid tangent from one constitutive trial evaluation per material
 * point.
 *
 * The supplied displacement field defines the trial current coordinates
 * directly as `x = X + u`; persistent model positions are therefore not needed
 * by this solver-facing evaluation. At every stiffness integration point,
 * Green-Lagrange strain is evaluated from the deformation gradient and the
 * constitutive law maps
 *
 *     (E_trial, state_committed) -> (S_trial, C_alg, state_trial).
 *
 * The resulting PK2 stress is an element-local intermediate used immediately for
 * the internal-force contribution
 *
 *     f_int = integral B^T S dV0.
 *
 * When `buffer != nullptr`, the same stress and tangent additionally form
 *
 *     K_mat = integral B^T C_alg B dV0,
 *
 * and the Total-Lagrangian geometric tangent assembled from
 * `grad(N_a)^T S grad(N_b)`. When `buffer == nullptr`, tangent assembly is
 * skipped after the same physical material update and internal-force assembly.
 *
 * @param buffer Optional dense tangent storage; null requests residual only.
 * @param nodal_forces Global nodal internal-force field to increment.
 * @param displacement Trial global nodal displacement field.
 * @return Mapped complete tangent, or an empty map for residual-only evaluation.
 */
template<Index N>
MapMatrix SolidElement<N>::stiffness_tangent(Precision*   buffer,
                                             NodeData&    nodal_forces,
                                             const Field& displacement) {
    logging::error(nodal_forces.components >= D,
        "SolidElement: nonlinear internal force requires at least three nodal components");

    const StaticMatrix<N, D> reference_coords   = this->node_coords_reference();
    const StaticMatrix<N, D> local_displacement = this->nodal_data<D>(displacement);
    const StaticMatrix<N, D> current_coords     = reference_coords + local_displacement;
    const auto& scheme                          = this->integration_scheme_stiffness();

    StaticMatrix<D * N, D * N> tangent = StaticMatrix<D * N, D * N>::Zero();

    // Every quadrature point performs exactly one physical constitutive update
    // from committed history into the corresponding persistent trial row.
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

        const Index      state_row = this->mp_index(ip);
        const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
        Precision*       new_state = &(*this->_model_data->material_state_new)(state_row, 0);

        VolumeStressPK2 stress;
        Mat6            material_tangent;
        evaluate_material(point.r, point.s, point.t, strain,
                          old_state, new_state, stress, material_tangent);

        const Precision measure = point.w * det0;

        // Internal force always comes from the same freshly evaluated PK2 stress.
        // No global integration-point stress scratch field is created or updated.
        const StaticVector<D * N> local_force = B.transpose() * stress.voigt() * measure;
        for (Index a = 0; a < N; ++a) {
            const Index node_id = static_cast<Index>(node_ids[a]);
            for (Dim d = 0; d < D; ++d) {
                nodal_forces(node_id, d) += local_force(D * a + d);
            }
        }

        // Residual-only evaluations still perform the complete physical material
        // trial update above, but skip all matrix assembly work.
        if (buffer == nullptr) {
            continue;
        }

        // Material contribution from the consistent algorithmic constitutive
        // tangent evaluated in the same trial state as the internal force.
        tangent.noalias() += measure * B.transpose() * material_tangent * B;

        // Total-Lagrangian geometric contribution generated by the local PK2
        // stress. Each scalar node-pair coefficient multiplies the 3x3 identity.
        const Mat3 S = stress.tensor();
        for (Index a = 0; a < N; ++a) {
            const Vec3 dNa = dN_dX.row(a).transpose();

            for (Index b = 0; b < N; ++b) {
                const Vec3 dNb = dN_dX.row(b).transpose();
                const Precision stress_coefficient = dNa.dot(S * dNb) * measure;

                for (Dim d = 0; d < D; ++d) {
                    tangent(D * a + d, D * b + d) += stress_coefficient;
                }
            }
        }
    }

    if (buffer == nullptr) {
        return MapMatrix(nullptr, 0, 0);
    }

    // Material and geometric tangents are analytically symmetric; remove only
    // floating-point asymmetry accumulated during quadrature and block assembly.
    tangent = Precision(0.5) * (tangent + tangent.transpose());

    MapMatrix mapped{buffer, D * N, D * N};
    mapped = tangent;
    return mapped;
}

/**
 * Recovers the global conductive heat flux at every solid integration point.
 *
 * The scalar nodal temperature field is differentiated with respect to global
 * reference coordinates. For the currently supported isotropic conductivity,
 * Fourier's law is
 *
 *     q = -k grad_X(T)
 *       = -k (dN/dX)^T T_e.
 *
 * Results are written to the element's globally enumerated `ELEMENT_IP` rows
 * using the same quadrature rule that defines `num_ip()`. The output components
 * are the global x-, y- and z-components of heat flux.
 *
 * @param heat_flux Global integration-point field receiving three heat-flux
 *                  components per row.
 * @param temperature Scalar global nodal temperature field.
 */
template<Index N>
void SolidElement<N>::compute_heat_flux(Field& heat_flux, const Field& temperature) {
    // Validate the input and output field layouts before gathering element data.
    logging::error(temperature.domain == FieldDomain::NODE,
        "SolidElement: temperature field must use the NODE domain");
    logging::error(temperature.components == 1,
        "SolidElement: temperature field must have exactly one component");
    logging::error(heat_flux.domain == FieldDomain::ELEMENT_IP,
        "SolidElement: heat-flux field must use the ELEMENT_IP domain");
    logging::error(heat_flux.components >= D,
        "SolidElement: heat-flux field requires at least three components");

    const auto& scheme = this->integration_scheme_stiffness();
    logging::error(this->ip_index(scheme.count()) <= heat_flux.rows,
        "SolidElement: heat-flux field is too small for element ", this->elem_id);

    // Gather the reference geometry, nodal temperatures and isotropic thermal
    // conductivity once for all integration points.
    const StaticMatrix<N, D> reference_coords  = this->node_coords_reference();
    const StaticVector<N>    local_temperature = this->template nodal_data<1>(temperature);
    auto*                    section           = this->get_section();

    logging::error(section->material_->has_thermal_conductivity(),
        "Material has no thermal conductivity at element ", this->elem_id);
    const Precision conductivity = section->material_->get_thermal_conductivity();

    // Apply Fourier's law at every globally enumerated thermal output point.
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        Precision det0 = Precision(0);
        const StaticMatrix<N, D> dN_dX = this->shape_derivatives_reference(
            reference_coords,
            point.r,
            point.s,
            point.t,
            det0
        );
        const Vec3 flux = -conductivity * dN_dX.transpose() * local_temperature;
        const Index row = this->ip_index(ip);

        for (Dim component = 0; component < D; ++component) {
            heat_flux(row, component) = flux(component);
        }
    }
}

/**
 * Computes the linear element compliance contribution `u^T K u`.
 *
 * The mechanical stiffness queried here is the reference-configuration linear
 * operator, so compliance evaluation does not advance persistent material state.
 *
 * @param displacement Global nodal displacement field.
 * @param result Element result field receiving the compliance contribution.
 */
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
 * The derivative is evaluated on the same reference-configuration small-strain
 * kinematics as `stiffness()`. Constitutive history is read from the committed
 * state and no persistent trial target is supplied.
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

    // Gather the element displacement and reference geometry once. Compliance
    // sensitivities follow the same linear small-strain kinematics as stiffness().
    auto local_disp_mat     = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    const auto reference_coords = this->node_coords_reference();
    const auto& scheme          = this->integration_scheme_stiffness();

    const Precision scaling    = element_stiffness_scale();
    Vec3            derivative = Vec3::Zero();

    // Integrate the energy sensitivity for all three rotation parameters.
    for (Index n = 0; n < scheme.count(); ++n) {
        const auto point = scheme.get_point(n);
        Precision det;

        const StaticMatrix<n_strain, D * N> B =
            this->strain_displacements(reference_coords, point.r, point.s, point.t, det);
        const StaticVector<n_strain> strain = B * local_displacement;

        // Evaluate transformation derivatives at the physical reference point
        // from committed material history without producing a trial state.
        const Vec3 position_reference =
            this->interpolate<D>(reference_coords, point.r, point.s, point.t);
        const Index      state_row = this->mp_index(n);
        const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);

        const auto tangent_derivatives = get_section()->tangent_rotation_derivatives(
            position_reference,
            additional_rotation,
            additional_rotation_derivatives,
            old_state,
            nullptr
        );

        // Accumulate epsilon^T C'_tan epsilon with the reference volume measure.
        for (Index i = 0; i < 3; ++i) {
            derivative(i) += scaling * point.w * strain.dot(tangent_derivatives[i] * strain) * det;
        }
    }

    result(elem_id, 0) = derivative(0);
    result(elem_id, 1) = derivative(1);
    result(elem_id, 2) = derivative(2);
}

} // namespace fem::model
