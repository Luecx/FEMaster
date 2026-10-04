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
 * Exact recovery evaluates F and PK2 stress at displacement. Affine recovery
 * evaluates them at linearization (nullptr means u_L = 0) and differentiates
 * both Green-Lagrange strain and the complete PK2-to-Cauchy push-forward.
 * Reference thermal strain contributes -C0 epsilon_th without advancing history.
 * Integration-point values are copied or extrapolated to topology node locations.
 *
 * @param strain Optional total strain output.
 * @param stress Optional Cauchy stress output.
 * @param displacement Requested nodal displacement.
 * @param rst Natural material-point or nodal output coordinates.
 * @param offset First output row belonging to this element.
 * @param thermal_free_strain Optional reference-only isotropic free strain.
 * @param linearization Affine expansion point; nullptr denotes zero displacement.
 */
template<Index N>
void SolidElement<N>::compute_stress_strain(
    Field*           strain,
    Field*           stress,
    const Field&     displacement,
    const RowMatrix& rst,
    int              offset,
    const Field*     linearization,
    const Field*     thermal_free_strain
) {
    const bool exact_state = linearization == &displacement;

    // Validate output coordinates and the supported thermal/kinematic combination
    logging::error(strain != nullptr || stress != nullptr,
        "SolidElement: compute_stress_strain requires at least one output field");
    logging::error(rst.cols() >= 3,
        "SolidElement: stress/strain evaluation coordinates require at least 3 columns");
    logging::error((!exact_state && linearization == nullptr) || thermal_free_strain == nullptr,
        "SolidElement: thermal free strain recovery is supported only for linear kinematics");
    logging::error(!thermal_free_strain || (thermal_free_strain->domain == FieldDomain::ELEMENT_NODAL && thermal_free_strain->components == 1),
        "SolidElement: thermal free strain must be scalar ELEMENT_NODAL data");

    const auto&     scheme       = this->integration_scheme_stiffness();
    const RowMatrix ip_rst       = this->stress_strain_ip_rst();
    const bool      output_at_ip = rst.rows() == ip_rst.rows() && rst.leftCols(3).isApprox(ip_rst);
    const bool output_at_nodes   = rst.rows() == static_cast<Eigen::Index>(N);

    logging::error(output_at_ip || output_at_nodes,
        "SolidElement: stress/strain output must use integration points or element nodes");

    const auto reference_coords       = this->node_coords_reference();
    const auto local_displacement     = this->nodal_data<3>(displacement);
    StaticMatrix<N, D> local_state = StaticMatrix<N, D>::Zero();
    if (exact_state) {
        local_state = local_displacement;
    } else if (linearization) {
        local_state = this->nodal_data<D>(*linearization);
    }
    const StaticMatrix<N, D> local_delta    = local_displacement - local_state;
    const StaticMatrix<N, D> current_coords = reference_coords + local_state;

    StaticVector<N> nodal_thermal_strain = StaticVector<N>::Zero();
    if (thermal_free_strain) {
        for (Index node = 0; node < N; ++node) {
            nodal_thermal_strain(node) = (*thermal_free_strain)(static_cast<Index>(this->elem_nodal_offset) + node, 0);
        }
    }

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
        const Mat3 F = this->deformation_gradient(reference_coords, current_coords, point.r, point.s, point.t);
        const VolumeStrainGreenLagrange green = VolumeStrainGreenLagrange::from_deformation_gradient(F);
        VolumeStressPK2 second_pk;
        Mat6            material_tangent;
        evaluate_material(
            point.r, point.s, point.t, green, old_state, nullptr, second_pk,
            exact_state ? nullptr : &material_tangent);

        Vec6 recovered_strain = green.voigt();
        const Mat3 sigma = second_pk.to_cauchy(F).tensor();
        Mat3 recovered_stress = sigma;
        if (!exact_state) {
            // Differentiate E and sigma = F S F^T / det(F) in the same direction
            const Mat3 delta_F = local_delta.transpose() * dN_dX;
            const Mat3 delta_E = Precision(0.5) * (F.transpose() * delta_F + delta_F.transpose() * F);
            const Vec6 delta_strain = VolumeStrainGreenLagrange(delta_E).voigt();
            recovered_strain += delta_strain;
            Vec6 mechanical_increment = delta_strain;
            if (thermal_free_strain) {
                const Precision free = this->shape_function(point.r, point.s, point.t).dot(nodal_thermal_strain);
                mechanical_increment.head<3>().array() -= free;
            }
            const Mat3 delta_S = VolumeStressPK2(Vec6(material_tangent * mechanical_increment)).tensor();
            const Mat3 S = second_pk.tensor();
            recovered_stress += (delta_F * S * F.transpose() + F * delta_S * F.transpose()
                               + F * S * delta_F.transpose()) / F.determinant()
                              - (F.inverse() * delta_F).trace() * sigma;
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
 * @param offset First element-nodal row belonging to this element.
 * @return True when the element uses J2 plasticity and contributes PEEQ.
 */
template<Index N>
bool SolidElement<N>::compute_peeq(Field& peeq, int offset) {
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
 * A null linearization denotes u_L = 0. Every request evaluates the same
 * Green-Lagrange/PK2 material law and its complete tangent K_T = K_M + K_G.
 * Forces away from u_L are the affine approximation f_L + K_T (u - u_L).
 * Trial history is written only for an exact physical evaluation at u_L.
 *
 * Reference thermal loading is an affine correction -B0^T C0 epsilon_th.
 * Reference geometric requests contract the linearly recovered PK2 prestress
 * with the reference gradients, preserving the initial-stress buckling operator.
 * Finite-state geometric requests use PK2 stress at the supplied u_L.
 *
 * @param tangent_buffer Optional storage for the complete element tangent.
 * @param geometric_tangent_buffer Optional storage for the geometric operator.
 * @param internal_force Optional global nodal force accumulator.
 * @param displacement State whose exact or affine internal force is requested.
 * @param linearization Tangent state; nullptr denotes the reference zero state.
 * @param thermal_free_strain Optional reference-only scalar element-nodal free strain.
 * @param update_state Write trial history for an exact evaluation at u_L.
 */
template<Index N>
MapMatrix SolidElement<N>::evaluate(
    Precision*   tangent_buffer,
    Precision*   geometric_tangent_buffer,
    NodeData*    internal_force,
    const Field* displacement,
    const Field* linearization,
    const Field* thermal_free_strain,
    bool         update_state
) {
    const bool with_tangent   = tangent_buffer != nullptr;
    const bool with_geometric = geometric_tangent_buffer != nullptr;
    const bool with_force     = internal_force != nullptr;
    if (!with_tangent && !with_geometric && !with_force) {
        return MapMatrix(nullptr, 0, 0);
    }

    // Validate force storage and keep affine queries state-neutral
    logging::error(!with_force || displacement != nullptr,
        "SolidElement: internal force evaluation requires displacement");
    logging::error(!update_state || (linearization != nullptr && displacement == linearization),
        "SolidElement: material state requires an exact evaluation at the linearization state");
    if (with_force) {
        logging::error(internal_force->components >= D,
            "SolidElement: internal force requires at least three nodal components");
    }
    logging::error(!thermal_free_strain || (linearization == nullptr && !update_state),
        "SolidElement: thermal free strain is supported only for reference linearization");
    logging::error(!thermal_free_strain || (thermal_free_strain->domain == FieldDomain::ELEMENT_NODAL && thermal_free_strain->components == 1),
        "SolidElement: thermal free strain must be scalar ELEMENT_NODAL data");

    // Gather the tangent state and the requested affine displacement increment
    const StaticMatrix<N, D> reference_coords = this->node_coords_reference();
    StaticMatrix<N, D> local_linearization = StaticMatrix<N, D>::Zero();
    StaticMatrix<N, D> local_delta         = StaticMatrix<N, D>::Zero();
    if (linearization) local_linearization = this->nodal_data<D>(*linearization);
    if (displacement)  local_delta = this->nodal_data<D>(*displacement) - local_linearization;
    const StaticMatrix<N, D> current_coords = reference_coords + local_linearization;
    const StaticMatrix<D, N> delta_matrix(local_delta.transpose());
    const auto delta = Eigen::Map<const StaticVector<D * N>>(delta_matrix.data(), D * N);

    StaticVector<N> nodal_thermal_strain = StaticVector<N>::Zero();
    if (thermal_free_strain) {
        for (Index node = 0; node < N; ++node) {
            nodal_thermal_strain(node) = (*thermal_free_strain)(this->elem_nodal_offset + node, 0);
        }
    }

    // Request the material tangent only for matrices, affine forces or thermal stress
    const bool affine_force = with_force && displacement != linearization;
    const bool need_tangent = with_tangent || affine_force;
    const bool need_material = need_tangent || thermal_free_strain
                            || (with_geometric && linearization == nullptr);
    StaticMatrix<D * N, D * N> tangent   = StaticMatrix<D * N, D * N>::Zero();
    StaticMatrix<D * N, D * N> geometric = StaticMatrix<D * N, D * N>::Zero();
    StaticVector<D * N> force           = StaticVector<D * N>::Zero();

    // Integrate the single finite-strain constitutive response over reference volume
    const auto& scheme = this->integration_scheme_stiffness();
    for (Index ip = 0; ip < scheme.count(); ++ip) {
        const auto point = scheme.get_point(ip);
        Precision det0;
        const StaticMatrix<N, D> dN_dX = this->shape_derivatives_reference(
            reference_coords, point.r, point.s, point.t, det0);
        const Mat3 F = this->deformation_gradient(reference_coords, current_coords, point.r, point.s, point.t);
        const StaticMatrix<n_strain, D * N> B = this->green_lagrange_strain_displacement(dN_dX, F);
        const VolumeStrainGreenLagrange strain = VolumeStrainGreenLagrange::from_deformation_gradient(F);

        const Index      state_row = this->mp_index(ip);
        const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
        Precision* new_state = update_state ? &(*this->_model_data->material_state_new)(state_row, 0) : nullptr;
        VolumeStressPK2 stress;
        Mat6            material_tangent;
        evaluate_material(
            point.r, point.s, point.t, strain, old_state, new_state, stress,
            need_material ? &material_tangent : nullptr);
        const Precision measure = point.w * det0;

        // Keep reference thermal strain as a linear correction to force and prestress
        Vec6 thermal_stress = Vec6::Zero();
        if (thermal_free_strain) {
            const Precision free = this->shape_function(point.r, point.s, point.t).dot(nodal_thermal_strain);
            thermal_stress = material_tangent * Vec6(free, free, free, 0, 0, 0);
        }
        if (with_force) {
            force.noalias() += measure * B.transpose() * (stress.voigt() - thermal_stress);
        }
        if (need_tangent) {
            tangent.noalias() += measure * B.transpose() * material_tangent * B;
        }

        // Reference buckling uses prestress linearized about zero, not a finite query at u
        const Mat3 S = stress.tensor();
        Mat3 geometric_stress = S;
        if (with_geometric && linearization == nullptr) {
            const Vec6 prestress = stress.voigt() + material_tangent * (B * delta) - thermal_stress;
            geometric_stress = VolumeStressPK2(prestress).tensor();
        }
        if (need_tangent || with_geometric) {
            for (Index a = 0; a < N; ++a) {
                const Vec3 dNa = dN_dX.row(a).transpose();
                for (Index b = 0; b < N; ++b) {
                    const Vec3 dNb = dN_dX.row(b).transpose();
                    const Precision full_coefficient = measure * dNa.dot(S * dNb);
                    const Precision geom_coefficient = measure * dNa.dot(geometric_stress * dNb);
                    for (Dim d = 0; d < D; ++d) {
                        if (need_tangent)   tangent(D * a + d, D * b + d) += full_coefficient;
                        if (with_geometric) geometric(D * a + d, D * b + d) += geom_coefficient;
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
    auto K = evaluate(buffer, nullptr, nullptr, nullptr, nullptr, nullptr, false);

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
