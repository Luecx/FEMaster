/**
 * @file element_solid_compute.ipp
 * @brief Implements solid stress recovery and nonlinear force/tangent assembly.
 *
 * Constitutive evaluations use globally enumerated committed/trial material-state
 * rows associated with each solid integration point. Result recovery and other
 * auxiliary evaluations are state-neutral. The physical nonlinear path always
 * evaluates PK2 stress and trial history; the material tangent is requested only
 * when the caller actually assembles a tangent matrix.
 *
 * @author Finn Eggers
 * @date 07.08.2026
 */

#pragma once

#include "../../cos/rectangular_system.h"
#include "../../material/isotropic_j2_elasticity.h"
#include "../../section/section_solid.h"

namespace fem::model {

/**
 * Evaluates Total-Lagrangian solid constitutive response with an optional tangent.
 *
 * The selected old/new state pointers are forwarded directly to the section.
 * PK2 stress is always evaluated and scaled by the element topology factor.
 * A null tangent is propagated through the section/material stack and therefore
 * avoids constitutive tangent work for
 * residual-only and nonlinear recovery evaluations.
 *
 * @param r First natural coordinate.
 * @param s Second natural coordinate.
 * @param t Third natural coordinate.
 * @param global_strain Green-Lagrange strain in global reference coordinates.
 * @param old_state Immutable material-point input state row.
 * @param new_state Optional material-point output state row.
 * @param global_stress PK2 stress returned in global reference coordinates.
 * @param global_tangent Optional consistent global material tangent `dS/dE`.
 */
template<Index N>
void SolidElement<N>::evaluate_material(
    Precision                        r,
    Precision                        s,
    Precision                        t,
    const VolumeStrainGreenLagrange&  global_strain,
    const Precision*                 old_state,
    Precision*                       new_state,
    VolumeStressPK2&                 global_stress,
    Mat6*                            global_tangent) {
    // Evaluate PK2 response in the global reference basis with optional tangent output
    get_section()->evaluate(
        material_position_reference(r, s, t),
        additional_material_rotation(),
        global_strain,
        old_state,
        new_state,
        global_stress,
        global_tangent
    );

    // Apply the element stiffness factor without changing constitutive history
    const Precision scaling = element_stiffness_scale();
    global_stress.voigt()   *= scaling;
    if (global_tangent != nullptr) {
        *global_tangent *= scaling;
    }
}

/**
 * Recovers exact or affine strain and Cauchy stress from the finite material law.
 *
 * Recovery uses the explicit target state (u,T) and base state (u0,T0).
 * Temperature is evaluated exactly at the fixed base geometry, while the
 * displacement response is linearized at the base state. At each integration
 * point this gives
 *
 *     S_target,0 = S(E(u0) - E_th(T)),
 *
 *     Delta S_u = C(u0,T0) Delta E.
 *
 * The Cauchy-stress continuation uses the target-temperature stress as its
 * anchor and the displacement derivative of the base state.
 *
 * Integration-point values are copied or extrapolated to topology node locations.
 *
 * @param strain Optional total strain output.
 * @param stress Optional Cauchy stress output.
 * @param target_displacement Requested nodal displacement state u.
 * @param target_temperature Requested nodal temperature state T.
 * @param rst Natural material-point or nodal output coordinates.
 * @param offset First output row belonging to this element.
 * @param base_displacement Base displacement state u0; nullptr denotes zero.
 * @param base_temperature Base temperature state T0; nullptr denotes the
 *        material stress-free temperature.
 */
template<Index N>
void SolidElement<N>::compute_stress_strain(
    Field*           strain,
    Field*           stress,
    const Field&     target_displacement,
    const Field*     target_temperature,
    const RowMatrix& rst,
    const Field*     base_displacement,
    const Field*     base_temperature
) {
    // First compiled output row belonging to this element.
    Index offset = static_cast<Index>(this->elem_nodal_offset);
    if ((strain && strain->domain == FieldDomain::ELEMENT_IP) || (stress && stress->domain == FieldDomain::ELEMENT_IP)) {
        offset = static_cast<Index>(this->elem_ip_offset);
    }

    const bool exact_displacement = base_displacement == &target_displacement;

    // Validate output coordinates
    logging::error(strain != nullptr || stress != nullptr,
        "SolidElement: compute_stress_strain requires at least one output field");
    logging::error(rst.cols() >= 3,
        "SolidElement: stress/strain evaluation coordinates require at least 3 columns");

    const auto&     scheme       = this->integration_scheme_stiffness();
    const RowMatrix ip_rst       = this->stress_strain_ip_rst();
    const bool      output_at_ip = rst.rows() == ip_rst.rows() && rst.leftCols(3).isApprox(ip_rst);
    const bool output_at_nodes   = rst.rows() == static_cast<Eigen::Index>(N);

    logging::error(output_at_ip || output_at_nodes,
        "SolidElement: stress/strain output must use integration points or element nodes");

    const auto reference_coords           = this->node_coords_reference();
    const auto local_target_displacement  = this->nodal_data<3>(target_displacement);
    StaticMatrix<N, D> local_state = StaticMatrix<N, D>::Zero();
    if (base_displacement) {
        local_state = this->nodal_data<D>(*base_displacement);
    }
    const StaticMatrix<N, D> local_delta    = local_target_displacement - local_state;
    const StaticMatrix<N, D> current_coords = reference_coords + local_state;

    // Build element-local free thermal strains for the base and target
    // temperature states. A null field denotes the material stress-free
    // temperature.
    auto material = this->material();
    const auto nodal_thermal_strain = [&](const Field* temperature_field) {
        StaticVector<N> result = StaticVector<N>::Zero();

        if (!temperature_field || !material || !material->has_thermal_expansion()) {
            return result;
        }

        const Precision T_zero = material->get_thermal_zero_temperature();
        const Precision alpha  = material->get_thermal_expansion();

        for (Index node = 0; node < N; ++node) {
            const Precision temperature =
                (*temperature_field)(static_cast<Index>(this->node_ids[node]), 0);
            result(node) = alpha * (temperature - T_zero);
        }

        return result;
    };

    const StaticVector<N> thermal_base   = nodal_thermal_strain(base_temperature);
    const StaticVector<N> thermal_target = nodal_thermal_strain(target_temperature);

    RowMatrix ip_strain = RowMatrix::Zero(scheme.count(), n_strain);
    RowMatrix ip_stress = RowMatrix::Zero(scheme.count(), n_strain);

    // Evaluate kinematics and material response only where constitutive state
    // actually exists: at the element's stiffness integration points.
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        const Index      state_row = this->mp_index(ip);
        const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);

        // Recover finite kinematics at the exact state or at the affine expansion point
        Precision det0;
        const StaticMatrix<N, D> dN_dX = this->shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const Mat3 F = this->deformation_gradient(
            reference_coords, current_coords, point.r, point.s, point.t);
        const VolumeStrainGreenLagrange green =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);

        const auto shape = this->shape_function(point.r, point.s, point.t);
        const Precision free_base   = shape.dot(thermal_base);
        const Precision free_target = shape.dot(thermal_target);

        // Tangent quantities are evaluated at (u0,T0).
        Vec6 constitutive_strain_base = green.voigt();
        constitutive_strain_base.head<3>().array() -= free_base;

        VolumeStressPK2 stress_base;
        Mat6            material_tangent;
        evaluate_material(
            point.r, point.s, point.t,
            VolumeStrainGreenLagrange(constitutive_strain_base),
            old_state, nullptr, stress_base,
            exact_displacement ? nullptr : &material_tangent);

        // The stress anchor for recovery is evaluated exactly at (u0,T).
        VolumeStressPK2 stress_target_base = stress_base;
        if (target_temperature != base_temperature) {
            Vec6 constitutive_strain_target = green.voigt();
            constitutive_strain_target.head<3>().array() -= free_target;

            evaluate_material(
                point.r, point.s, point.t,
                VolumeStrainGreenLagrange(constitutive_strain_target),
                old_state, nullptr, stress_target_base, nullptr);
        }

        Vec6 recovered_strain = green.voigt();
        const Mat3 sigma_base   = stress_base.to_cauchy(F).tensor();
        const Mat3 sigma_target = stress_target_base.to_cauchy(F).tensor();
        Mat3 recovered_stress   = sigma_target;
        if (!exact_displacement) {
            // Differentiate E and sigma = F S F^T / det(F) in the same direction
            const Mat3 delta_F = local_delta.transpose() * dN_dX;
            const Mat3 delta_E = Precision(0.5) * (F.transpose() * delta_F + delta_F.transpose() * F);
            const Vec6 delta_strain = VolumeStrainGreenLagrange(delta_E).voigt();
            recovered_strain += delta_strain;
            // Temperature is fixed during the displacement perturbation:
            //
            //     Delta E_th   = 0,
            //     Delta E_mech = Delta E.
            const Vec6 mechanical_increment = delta_strain;
            const Mat3 delta_S =
                VolumeStressPK2(Vec6(material_tangent * mechanical_increment)).tensor();
            const Mat3 S = stress_base.tensor();
            recovered_stress += (delta_F * S * F.transpose() + F * delta_S * F.transpose()
                               + F * S * delta_F.transpose()) / F.determinant()
                              - (F.inverse() * delta_F).trace() * sigma_base;
        }
        ip_strain.row(ip) = recovered_strain.transpose();
        ip_stress.row(ip) = VolumeStressCauchy(recovered_stress).voigt().transpose();
    }

    // Integration-point output is the constitutive result itself and must not be
    // projected through a recovery basis.
    if (output_at_ip) {
        for (Eigen::Index n = 0; n < rst.rows(); ++n) {
            const Index row = static_cast<Index>(offset + n);
            for (Dim component = 0; component < n_strain; ++component) {
                if (strain) (*strain)(row, component) = ip_strain(n, component);
                if (stress) (*stress)(row, component) = ip_stress(n, component);
            }
        }
        return;
    }

    // Nodal values are reconstructed from the integration-point samples in
    // natural coordinates using the topology-specific constant operator.
    const RowMatrix& E      = this->extrapolation_matrix();
    logging::error(E.rows() == rst.rows() && E.cols() == static_cast<Eigen::Index>(scheme.count()),
        "SolidElement: invalid extrapolation matrix for element ", this->elem_id);

    const RowMatrix nodal_strain = E * ip_strain;
    const RowMatrix nodal_stress = E * ip_stress;

    for (Eigen::Index n = 0; n < rst.rows(); ++n) {
        const Index row = static_cast<Index>(offset + n);
        for (Dim component = 0; component < n_strain; ++component) {
            if (strain) (*strain)(row, component) = nodal_strain(n, component);
            if (stress) (*stress)(row, component) = nodal_stress(n, component);
        }
    }
}

/**
 * Recovers accumulated equivalent plastic strain from committed solid material
 * history and extrapolates the integration-point values to the element nodes.
 *
 * Each solid integration point owns one material-state row. J2 materials expose
 * PEEQ directly from that committed row, so recovery does not reevaluate the
 * constitutive law or modify trial history. The same topology-specific
 * extrapolation matrix used by stress recovery maps the scalar IP values to
 * element-nodal output.
 *
 * A J2 material whose nonlinear state storage has not been initialized yet
 * contributes zero PEEQ at all element nodes. Non-J2 materials return false and
 * are excluded from the subsequent model-wide nodal average.
 *
 * @param peeq Scalar ELEMENT_NODAL output field.
 * @return True when the element uses J2 plasticity and contributes PEEQ.
 */
template<Index N>
bool SolidElement<N>::compute_peeq(Field& peeq) {
    // First element-nodal result row belonging to this element.
    Index offset = static_cast<Index>(this->elem_nodal_offset);

    logging::error(peeq.domain == FieldDomain::ELEMENT_NODAL && peeq.components == 1,
        "SolidElement: PEEQ recovery requires scalar ELEMENT_NODAL output");

    auto mat = get_section()->material_;
    if (!mat || !mat->has_elasticity()) return false;

    const auto* j2 = mat->elasticity()->template as<material::IsotropicJ2Elasticity>();
    if (!j2) return false;

    // A J2 material that has not entered a nonlinear analysis is still virgin.
    const auto& state = this->_model_data->material_state_old;
    if (!state || state->components < j2->state_size()) {
        for (Index node = 0; node < N; ++node) {
            peeq(static_cast<Index>(offset) + node, 0) = Precision(0);
        }
        return true;
    }

    // Read the accepted history at every constitutive integration point.
    const auto& scheme = this->integration_scheme_stiffness();
    RowMatrix ip_peeq  = RowMatrix::Zero(scheme.count(), 1);

    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const Index state_row      = this->mp_index(ip);
        const Precision* old_state = &(*state)(state_row, 0);
        ip_peeq(ip, 0)             = j2->equivalent_plastic_strain(old_state);
    }

    // Extrapolate the scalar integration-point history to the element nodes.
    const RowMatrix& E      = this->extrapolation_matrix();
    logging::error(E.rows() == N && E.cols() == static_cast<Eigen::Index>(scheme.count()),
        "SolidElement: invalid PEEQ extrapolation matrix for element ", this->elem_id);

    const RowMatrix nodal_peeq = E * ip_peeq;
    for (Index node = 0; node < N; ++node) {
        peeq(static_cast<Index>(offset) + node, 0) = nodal_peeq(node, 0);
    }

    return true;
}

/**
 * Evaluates the Total-Lagrangian element response at one linearization state.
 *
 * A null linearization denotes u0 = 0. Every request evaluates the same
 * Green-Lagrange/PK2 material law and the complete tangent K_T = K_M + K_G at
 * u0. Forces away from u0 are the affine approximation
 * f_int(u0) + K_T(u0) (u - u0). Trial history is written only for an exact
 * physical evaluation at u0.
 *
 * The separately requested geometric operator is different from the geometric
 * part already contained in K_T(u0): it is assembled only from the linearized
 * PK2 stress increment produced by the perturbation from u0 to u. Existing
 * stress at u0 therefore remains part of the complete tangent and is not repeated
 * in the separate geometric operator.
 *
 * @param tangent_buffer Optional storage for the complete element tangent at u0.
 * @param geometric_tangent_buffer Optional storage for the geometric stiffness
 *        generated by the linearized stress increment from u0 to u.
 * @param internal_force Optional global nodal force accumulator.
 * @param target_displacement Requested displacement state u.
 * @param target_temperature Requested temperature state T.
 * @param base_displacement Base displacement state u0; nullptr denotes zero.
 * @param base_temperature Base temperature state T0; nullptr denotes the
 *        material stress-free temperature.
 * @param update_state Write trial history for an exact evaluation at u_L.
 */
template<Index N>
MapMatrix SolidElement<N>::evaluate(
    Precision*   tangent_buffer,
    Precision*   geometric_tangent_buffer,
    NodeData*    internal_force,
    const Field* target_displacement,
    const Field* target_temperature,
    const Field* base_displacement,
    const Field* base_temperature,
    bool         update_state
) {
    // Requested element outputs:
    // - complete tangent stiffness matrix K_T(u0) = K_M(u0) + K_G(S0)
    // - geometric stiffness generated by the linearized stress increment from u0 to u
    // - internal force vector at u
    const bool with_tangent   = tangent_buffer           != nullptr;
    const bool with_geometric = geometric_tangent_buffer != nullptr;
    const bool with_force     = internal_force           != nullptr;

    // Validate requested outputs and evaluation state
    logging::error(with_tangent || with_geometric || with_force,
        "SolidElement: evaluation requires at least one requested output");
    logging::error(!with_force || target_displacement != nullptr,
        "SolidElement: internal force evaluation requires displacement");
    logging::error(!with_force || internal_force->components >= D,
        "SolidElement: internal force requires at least three nodal components");
    logging::error(!update_state ||
        (base_displacement != nullptr && target_displacement == base_displacement &&
         target_temperature == base_temperature),
        "SolidElement: material state requires an exact evaluation at the base state");

    // Gather the reference nodal coordinates and prepare the displacement states
    // used for the Total-Lagrangian linearization.
    StaticMatrix<N, D> reference_coords = this->node_coords_reference();
    StaticMatrix<N, D> disp_linear_base = StaticMatrix<N, D>::Zero();
    StaticMatrix<N, D> disp_delta       = StaticMatrix<N, D>::Zero();

    // Collect the displacement u0 at the linearization point and the displacement
    // increment du = u - u0.
    if (base_displacement)   disp_linear_base = this->nodal_data<D>(*base_displacement);
    if (target_displacement) disp_delta       = this->nodal_data<D>(*target_displacement) - disp_linear_base;

    // Build the nodal coordinates of the linearization configuration x0 = X + u0.
    StaticMatrix<N, D> linearization_coords = reference_coords + disp_linear_base;

    // Flatten the nodal displacement increment into the element DOF ordering.
    StaticMatrix<D, N> delta_matrix(disp_delta.transpose());
    const auto delta = Eigen::Map<const StaticVector<D * N>>(delta_matrix.data(), D * N);

    // Evaluate isotropic free thermal strain for both temperature states.
    // Null temperature fields denote the material stress-free temperature.
    auto material = this->material();
    const auto nodal_thermal_strain = [&](const Field* temperature_field) {
        StaticVector<N> result = StaticVector<N>::Zero();

        if (!temperature_field || !material || !material->has_thermal_expansion()) {
            return result;
        }

        const Precision T_zero = material->get_thermal_zero_temperature();
        const Precision alpha  = material->get_thermal_expansion();

        for (Index node = 0; node < N; ++node) {
            const Precision temperature =
                (*temperature_field)(static_cast<Index>(this->node_ids[node]), 0);
            result(node) = alpha * (temperature - T_zero);
        }

        return result;
    };

    const StaticVector<N> thermal_base   = nodal_thermal_strain(base_temperature);
    const StaticVector<N> thermal_target = nodal_thermal_strain(target_temperature);

    // -----------------------------------------------------------------------------
    // main loop
    // -----------------------------------------------------------------------------
    // Determine which quantities must be assembled internally.
    //
    // An affine force evaluation is required when the requested displacement u
    // differs from the linearization state u0. The internal force is then
    // approximated by
    //
    //     f_int(u) = f_int(u0) + K_T(u0) (u - u0).
    //
    // If u = u0, the internal force is evaluated directly at the linearization
    // state and no affine correction is required.
    const bool affine_force = with_force && target_displacement != base_displacement;

    // The complete tangent K_T = K_M + K_G is required either because it was
    // explicitly requested or because the affine force evaluation needs
    // K_T(u0) (u - u0).
    const bool need_tangent = with_tangent || affine_force;

    // The constitutive tangent C = dS/dE is required for:
    // - the material tangent contribution K_M = ∫ B^T C B dV,
    // - the stress increment used by the separately requested geometric stiffness.
    const bool need_material = need_tangent || with_geometric;

    // Initialize the element quantities accumulated over all integration points.
    StaticMatrix<D * N, D * N> tangent   = StaticMatrix<D * N, D * N>::Zero();
    StaticMatrix<D * N, D * N> geometric = StaticMatrix<D * N, D * N>::Zero();
    StaticVector<D * N>        force     = StaticVector<D * N>::Zero();

    // Integrate the single finite-strain constitutive response over reference volume
    const auto& scheme = this->integration_scheme_stiffness();
    for (Index ip = 0; ip < scheme.count(); ++ip) {

        // -----------------------------------------------------------------------------
        // strain-displacement / kinematics
        // -----------------------------------------------------------------------------
        const auto point = scheme.get_point(ip);
        // Transform the shape-function derivatives from natural coordinates (r,s,t) to
        // global reference coordinates (x,y,z) in undeformed (reference) space.
        Precision det0;
        const StaticMatrix<N, D> dN_dX = this->shape_derivatives_reference(reference_coords, point.r, point.s, point.t, det0);

        // Evaluate the deformation gradient F = dx0/dX at the linearization
        // configuration x0 = X + u0.
        const Mat3 F = this->deformation_gradient(reference_coords, linearization_coords, point.r, point.s, point.t);

        // Build the derivative of the Green-Lagrange strain with respect to the
        // element nodal displacements at the linearization state,
        //
        //     B = dE/du |u0.
        const StaticMatrix<n_strain, D * N> B = this->green_lagrange_strain_displacement(dN_dX, F);

        // Evaluate the Green-Lagrange strain at the linearization state,
        //
        //     E0 = 1/2 (F^T F - I).
        const VolumeStrainGreenLagrange strain =
            VolumeStrainGreenLagrange::from_deformation_gradient(F);

        const auto shape = this->shape_function(point.r, point.s, point.t);
        const Precision thermal_strain_base   = shape.dot(thermal_base);
        const Precision thermal_strain_target = shape.dot(thermal_target);

        // The complete tangent is evaluated at the base state (u0,T0).
        Vec6 constitutive_strain_base = strain.voigt();
        constitutive_strain_base.head<3>().array() -= thermal_strain_base;

        // -----------------------------------------------------------------------------
        // material evaluation
        // -----------------------------------------------------------------------------
        // Locate the material-point state belonging to the current integration point
        // and access its constitutive history. old_state contains the previously
        // committed material variables, while new_state provides writable storage
        // for the updated state when update_state is enabled.
        const Index      state_row = this->mp_index(ip);
        const Precision* old_state =                &(*this->_model_data->material_state_old)(state_row, 0);
        Precision*       new_state = update_state ? &(*this->_model_data->material_state_new)(state_row, 0) : nullptr;

        // Evaluate the constitutive response at the current integration point.
        // stress stores the second Piola-Kirchhoff stress S associated with the
        // mechanical strain E_mech,0. material_tangent optionally stores the
        // constitutive tangent C = dS/dE when required by a tangent operator.
        VolumeStressPK2 stress_base;
        Mat6            material_tangent;
        evaluate_material(
            point.r, point.s, point.t,
            VolumeStrainGreenLagrange(constitutive_strain_base),
            old_state, new_state, stress_base,
            need_material ? &material_tangent : nullptr);

        // The force anchor is evaluated at the target temperature but at the
        // same base geometry u0.
        VolumeStressPK2 stress_target_base = stress_base;
        if ((with_force || with_geometric) && target_temperature != base_temperature) {
            Vec6 constitutive_strain_target = strain.voigt();
            constitutive_strain_target.head<3>().array() -= thermal_strain_target;

            evaluate_material(
                point.r, point.s, point.t,
                VolumeStrainGreenLagrange(constitutive_strain_target),
                old_state, nullptr, stress_target_base, nullptr);
        }

        // -----------------------------------------------------------------------------
        // contributions to default stiffness and force vector
        // -----------------------------------------------------------------------------
        if (with_force)   force.noalias()   += det0 * point.w * B.transpose() * stress_target_base.voigt();
        if (need_tangent) tangent.noalias() += det0 * point.w * B.transpose() * material_tangent * B;

        // -----------------------------------------------------------------------------
        // geometric stiffness
        // -----------------------------------------------------------------------------
        //
        // In the Total-Lagrangian formulation, the geometric stiffness is obtained
        // from a PK2 stress state and the shape-function gradients in the reference
        // configuration:
        //
        //     K_G,ab = ∫ (grad_X N_a)^T S (grad_X N_b) I dV0.
        //
        // The geometric contribution to the complete tangent uses the actual stress S0
        // at the linearization state u0.
        //
        // A separately requested geometric stiffness uses only the stress increment
        // from u0 to the requested state u:
        //
        //     Delta S = C0 B0 (u - u0).
        //
        // Temperature is held fixed during this displacement perturbation. Thus,
        // the complete base stress S0, including thermal stress, contributes to
        // K_T(u0), while the separate geometric operator contains only Delta S.
        const Mat3 stress_linear_base = stress_base.tensor();

        Vec6 stress_increment = Vec6::Zero();
        if (with_geometric) {
            stress_increment =
                stress_target_base.voigt() - stress_base.voigt()
                + material_tangent * (B * delta);
        }

        const Mat3 stress_geometric = VolumeStressPK2(stress_increment).tensor();
        // Integrate the geometric contributions over the reference volume.
        if (need_tangent || with_geometric) {
            for (Index a = 0; a < N; ++a) {
                const Vec3 dNa = dN_dX.row(a).transpose();

                for (Index b = 0; b < N; ++b) {
                    const Vec3 dNb = dN_dX.row(b).transpose();

                    const Precision tangent_geometric_coefficient = det0 * point.w * dNa.dot(stress_linear_base * dNb);
                    const Precision geometric_coefficient         = det0 * point.w * dNa.dot(stress_geometric   * dNb);

                    for (Dim d = 0; d < D; ++d) {
                        if (need_tangent)   tangent  (D * a + d, D * b + d) += tangent_geometric_coefficient;
                        if (with_geometric) geometric(D * a + d, D * b + d) += geometric_coefficient;
                    }
                }
            }
        }
    }

    // Apply the affine increment only after assembling the complete tangent
    if (affine_force) force.noalias() += tangent * delta;
    if (with_force) {
        for (Index a = 0; a < N; ++a) {
            for (Dim d = 0; d < D; ++d) {
                (*internal_force)(node_ids[a], d) += force(D * a + d);
            }
        }
    }

    // Copy requested operators into caller-owned contiguous storage
    if (with_tangent) {
        MapMatrix mapped(tangent_buffer, D * N, D * N);
        mapped = tangent;
    }
    if (with_geometric) {
        MapMatrix mapped(geometric_tangent_buffer, D * N, D * N);
        mapped = geometric;
    }
    return with_tangent ? MapMatrix(tangent_buffer, D * N, D * N)
         : with_geometric ? MapMatrix(geometric_tangent_buffer, D * N, D * N)
                          : MapMatrix(nullptr, 0, 0);
}

/**
 * Recovers conductive heat flux as one vector per element node.
 *
 * Fourier's law is evaluated from the scalar nodal temperature interpolation,
 *
 * q = -k grad_X(T)
 * = -k (dN/dX)^T T_e.
 *
 * The gradient is first evaluated at the formulation's stable stiffness
 * integration points and then recovered to the element nodes through the same
 * topology-specific extrapolation operator used by structural nodal recovery.
 * This is important for reduced and degenerated elements, where direct
 * differentiation at a geometric node can be singular or intentionally differs
 * from the formulation's recovery basis.
 *
 * No global `ELEMENT_IP` heat-flux field is created. The temporary integration-
 * point values remain local to this function and the final vectors are written
 * directly into the element's disjoint `ELEMENT_NODAL` range. This makes the
 * operation safe for model-level parallel recovery.
 *
 * @param heat_flux Global element-nodal field receiving three heat-flux
 * components per element node.
 * @param temperature Scalar global nodal temperature field.
 */
template<Index N>
void SolidElement<N>::compute_heat_flux(Field& heat_flux, const Field& temperature) {
    // Validate the scalar primary field and element-nodal output layout before
    // gathering any element-local data.
    logging::error(temperature.domain == FieldDomain::NODE,
        "SolidElement: temperature field must use the NODE domain");
    logging::error(temperature.components == 1,
        "SolidElement: temperature field must have exactly one component");
    logging::error(heat_flux.domain == FieldDomain::ELEMENT_NODAL,
        "SolidElement: heat-flux field must use the ELEMENT_NODAL domain");
    logging::error(heat_flux.components >= D,
        "SolidElement: heat-flux field requires at least three components");

    const Index offset = static_cast<Index>(this->elem_nodal_offset);
    logging::error(offset + N <= heat_flux.rows,
        "SolidElement: heat-flux field is too small for element ", this->elem_id);

    // Gather reference geometry, scalar nodal temperatures and conductivity once
    // for the complete recovery operation.
    const StaticMatrix<N, D> reference_coords  = this->node_coords_reference();
    const StaticVector<N>    local_temperature = this->template nodal_data<1>(temperature);
    auto*                    section           = this->get_section();

    logging::error(section->material_->has_thermal_conductivity(),
        "Material has no thermal conductivity at element ", this->elem_id);
    const Precision conductivity = section->material_->get_thermal_conductivity();

    // Evaluate Fourier's law only in local temporary storage at the formulation's
    // stable integration points. Nothing is written to global ELEMENT_IP storage.
    const auto& scheme = this->integration_scheme_stiffness();
    RowMatrix ip_flux  = RowMatrix::Zero(scheme.count(), D);

    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);

        Precision det0                 = Precision(0);
        const StaticMatrix<N, D> dN_dX = this->shape_derivatives_reference(
            reference_coords,
            point.r,
            point.s,
            point.t,
            det0
        );

        const Vec3 flux = -conductivity * dN_dX.transpose() * local_temperature;
        for (Dim component = 0; component < D; ++component) {
            ip_flux(ip, component) = flux(component);
        }
    }

    // Recover the integration-point vectors to all element nodes. The constant
    // topology-specific matrix also handles reduced-integration formulations by
    // reproducing their admissible recovery space instead of evaluating a
    // potentially singular nodal Jacobian directly.
    const RowMatrix& E      = this->extrapolation_matrix();
    logging::error(E.rows() == static_cast<Eigen::Index>(N) && E.cols() == static_cast<Eigen::Index>(scheme.count()),
        "SolidElement: invalid heat-flux extrapolation matrix for element ", this->elem_id);

    const RowMatrix nodal_flux = E * ip_flux;
    logging::error(nodal_flux.allFinite(),
        "SolidElement: nodal heat flux contains NaN or Inf at element ", this->elem_id);

    // Write only the element-owned disjoint rows. Model-level OpenMP recovery can
    // therefore evaluate different elements without synchronization.
    for (Index local_node = 0; local_node < N; ++local_node) {
        const Index row = offset + local_node;
        for (Dim component = 0; component < D; ++component) {
            heat_flux(row, component) = nodal_flux(local_node, component);
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
    // Obtain the state-neutral reference tangent through the common mechanical path
    Precision buffer[D * N * D * N];
    auto K = evaluate(buffer, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, false);

    auto local_disp_mat     = StaticMatrix<3, N>(this->nodal_data<3>(displacement).transpose());
    auto local_displacement = Eigen::Map<StaticVector<3 * N>>(local_disp_mat.data(), 3 * N);

    Precision strain_energy = local_displacement.dot((K * local_displacement));
    result(elem_id, 0)      = strain_energy;
}

/**
 * Computes the compliance derivative with respect to the three additional
 * material-orientation angles from the derivative of the constitutive tangent.
 *
 * With equilibrium `K u = f` and compliance `J = f^T u`, differentiation gives
 *
 * J' = -u^T K' u
 * = -integral eps^T C'_tan eps dV.
 *
 * The derivative is evaluated on the same reference-configuration small-strain
 * kinematics as the reference branch of `evaluate()`. Constitutive history is
 * read from the committed state and no persistent trial target is supplied.
 *
 * @param displacement Current nodal displacement field.
 * @param result Element result field receiving the three angle derivatives.
 */
template<Index N>
void SolidElement<N>::compute_compliance_angle_derivative(Field& displacement, Field& result) {
    if (!this->_model_data || !this->_model_data->material_orientation) {
        return;
    }

    auto angles_field                       = this->_model_data->material_orientation;
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
    // sensitivities use the tangent of evaluate() at zero displacement.
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

        // Form the compatible small-strain field in the reference configuration
        const StaticMatrix<N, D> dN_dX = this->shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det);
        const StaticMatrix<n_strain, D * N> B = this->strain_displacement(dN_dX);
        const StaticVector<n_strain> strain   = B * local_displacement;

        // Evaluate transformation derivatives at the physical reference point
        // from committed material history without producing a trial state.
        const Vec3 position_reference = this->interpolate<D>(reference_coords, point.r, point.s, point.t);
        const Index      state_row    = this->mp_index(n);
        const Precision* old_state    = &(*this->_model_data->material_state_old)(state_row, 0);

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
