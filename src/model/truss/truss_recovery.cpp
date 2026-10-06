/**
 * @file truss_recovery.cpp
 * @brief Implements T3 stress, history, compliance and section-force recovery.
 *
 * Recovery is deliberately state-neutral. Constitutive calls read committed
 * material history but never write trial history. The requested displacement
 * field defines the state being recovered explicitly; no mechanical result is
 * inferred indirectly from the model POSITION field.
 *
 * @see T3
 *
 * @author Finn Eggers
 * @date 04.10.2026
 */

#include "truss.h"

#include "../../material/isotropic_j2_elasticity.h"

namespace fem {
namespace model {

RowMatrix T3::stress_strain_nodal_rst() {
    RowMatrix rst(N, 3);
    rst.setZero();

    // Natural coordinates of the two-node line interpolation:
    //
    //     node 1: r = -1,
    //     node 2: r = +1.
    rst(0, 0) = Precision(-1);
    rst(1, 0) = Precision( 1);
    return rst;
}

RowMatrix T3::stress_strain_ip_rst() {
    // The one-point axial state is located at the element midpoint r = 0.
    RowMatrix rst(1, 3);
    rst.setZero();
    return rst;
}

/**
 * Evaluates the axial strain and physical Cauchy stress at the requested state.
 *
 * The optional linearization field defines the exact base configuration u0.
 * Without one, the undeformed reference configuration is used,
 *
 *     u0 = 0.
 *
 * The exact configuration at u0 is described by
 *
 *     r0      = (X2 + u02) - (X1 + u01),
 *     lambda0 = ||r0|| / L0,
 *     n0      = r0 / ||r0||.
 *
 * The corresponding Green-Lagrange strain is
 *
 *     E0 = 1/2 (lambda0^2 - 1).
 *
 * Thermal expansion contributes an eigenstrain
 *
 *     E_th = alpha (T - T0),
 *
 * where the truss temperature is taken as the mean nodal temperature. The
 * constitutive strain at the linearization point is therefore
 *
 *     E_mech,0 = E0 - E_th.
 *
 * The constitutive model returns the work-conjugate PK2 stress S0 and material
 * tangent C0 at this exact base state.
 *
 * The displacement perturbation between the requested state u and the exact
 * base state u0 is
 *
 *     Delta u = u - u0,
 *
 * giving
 *
 *     Delta r = Delta u2 - Delta u1.
 *
 * Linearizing stretch and Green-Lagrange strain about u0 gives
 *
 *     Delta lambda = n0 . Delta r / L0,
 *
 *     Delta E      = lambda0 Delta lambda
 *                  = lambda0 n0 . Delta r / L0.
 *
 * Strain and PK2 stress are continued linearly from the exact base state,
 *
 *     E ~= E0 + Delta E,
 *
 *     S ~= S0 + C0 Delta E.
 *
 * Physical axial stress follows from
 *
 *     sigma = lambda S.
 *
 * Its first-order expansion about u0 is therefore
 *
 *     sigma ~= lambda0 S0
 *            + S0 Delta lambda
 *            + lambda0 Delta S.
 *
 * The product Delta lambda Delta S is second order and is deliberately
 * omitted.
 *
 * @param displacement Requested displacement state u.
 * @param linearization Optional exact base state u0; null denotes u0 = 0.
 * @return Pair containing axial strain and physical axial stress.
 */
std::pair<Precision, Precision> T3::evaluate_axial_response(
    const Field& target_displacement,
    const Field* target_temperature,
    const Field* base_displacement,
    const Field* base_temperature
) {
    const Vec3 X1 = node_position_reference(0);
    const Vec3 X2 = node_position_reference(1);
    const Vec3 u1 = target_displacement.row_vec3(static_cast<Index>(node_ids[0]));
    const Vec3 u2 = target_displacement.row_vec3(static_cast<Index>(node_ids[1]));

    Vec3 u01 = Vec3::Zero();
    Vec3 u02 = Vec3::Zero();

    if (base_displacement) {
        u01 = base_displacement->row_vec3(static_cast<Index>(node_ids[0]));
        u02 = base_displacement->row_vec3(static_cast<Index>(node_ids[1]));
    }

    const Vec3      reference_axis = X2 - X1;
    const Precision L0             = reference_axis.norm();
    const Vec3      axis_base      = (X2 + u02) - (X1 + u01);
    const Precision length_base    = axis_base.norm();

    auto elasticity = get_elasticity();

    logging::error(L0 > Precision(0),
        "T3: zero reference length in evaluate_axial_response for element ", this->elem_id);
    logging::error(length_base > Precision(0),
        "T3: zero length at the linearization state in evaluate_axial_response for element ",
        this->elem_id);

    // -------------------------------------------------------------------------
    // exact state at the linearization point u0
    // -------------------------------------------------------------------------

    // The exact configuration at u0 is described by
    //
    //     r0      = (X2 + u02) - (X1 + u01),
    //     lambda0 = ||r0|| / L0,
    //     n0      = r0 / ||r0||.
    //
    // The corresponding Green-Lagrange strain is
    //
    //     E0 = 1/2 (lambda0^2 - 1).
    const Precision lambda0 = length_base / L0;
    const Vec3      n0      = axis_base / length_base;

    const AxialStrainGreenLagrange strain_base =
        AxialStrainGreenLagrange::from_stretch(lambda0);

    // -------------------------------------------------------------------------
    // thermal contribution
    // -------------------------------------------------------------------------

    // Thermal expansion is evaluated for both temperature states. A null field
    // denotes the material stress-free temperature.
    const auto material = get_material();
    const auto thermal_strain = [&](const Field* temperature_field) {
        if (!temperature_field || !material->has_thermal_expansion()) {
            return Precision(0);
        }

        const Precision T_zero = material->get_thermal_zero_temperature();
        const Precision alpha  = material->get_thermal_expansion();

        Precision temperature = Precision(0);
        for (Index node = 0; node < N; ++node) {
            temperature += (*temperature_field)(
                static_cast<Index>(node_ids[node]), 0);
        }

        temperature /= static_cast<Precision>(N);
        return alpha * (temperature - T_zero);
    };

    const AxialStrainGreenLagrange mechanical_strain_base(
        strain_base.value() - thermal_strain(base_temperature));
    const AxialStrainGreenLagrange mechanical_strain_target(
        strain_base.value() - thermal_strain(target_temperature));

    // -------------------------------------------------------------------------
    // material model evaluation
    // -------------------------------------------------------------------------

    AxialStressPK2 stress_base;
    AxialStressPK2 stress_target_base;
    Precision      material_tangent = Precision(0);

    const Index      state_row = this->mp_index(0);
    const Precision* old_state =
        &(*this->_model_data->material_state_old)(state_row, 0);

    elasticity->evaluate(
        mechanical_strain_base,
        old_state,
        nullptr,
        stress_base,
        &material_tangent
    );

    if (target_temperature == base_temperature) {
        stress_target_base = stress_base;
    } else {
        elasticity->evaluate(
            mechanical_strain_target,
            old_state,
            nullptr,
            stress_target_base,
            nullptr
        );
    }

    // -------------------------------------------------------------------------
    // linear continuation from u0 to u
    // -------------------------------------------------------------------------

    // The displacement perturbation is
    //
    //     Delta u = u - u0
    //
    // and therefore
    //
    //     Delta r = Delta u2 - Delta u1.
    const Vec3 delta_axis = (u2 - u02) - (u1 - u01);

    // Linearizing the stretch and Green-Lagrange strain about u0 gives
    //
    //     Delta lambda = n0 . Delta r / L0,
    //
    //     Delta E      = lambda0 Delta lambda
    //                  = lambda0 n0 . Delta r / L0.
    const Precision delta_lambda = n0.dot(delta_axis) / L0;
    const Precision delta_strain = lambda0 * delta_lambda;

    // Continue strain and PK2 stress linearly from the exact base state:
    //
    //     E ~= E0 + Delta E,
    //
    //     S ~= S0 + C0 Delta E.
    const Precision strain_value     = strain_base.value() + delta_strain;
    const Precision stress_increment = material_tangent * delta_strain;

    // Physical axial stress is
    //
    //     sigma = lambda S.
    //
    // Its first-order expansion about u0 is therefore
    //
    //     sigma ~= lambda0 S0
    //            + S0 Delta lambda
    //            + lambda0 Delta S.
    //
    // The product Delta lambda Delta S is second order and is deliberately
    // omitted.
    const Precision stress_value = lambda0 * stress_target_base.value()
                                 + stress_base.value() * delta_lambda
                                 + lambda0 * stress_increment;

    return {strain_value, stress_value};
}

/**
 * Recovers axial strain and physical Cauchy stress.
 *
 * Reference recovery uses the infinitesimal axial strain
 *
 *     epsilon = N0 . (u2 - u1) / L0,
 *
 * where
 *
 *     N0 = (X2 - X1) / L0.
 *
 * Exact finite recovery constructs the deformed axis directly from the supplied
 * displacement:
 *
 *     r = (X2 + u2) - (X1 + u1),
 *     l = ||r||,
 *     lambda = l / L0,
 *
 * and evaluates the Total-Lagrangian strain
 *
 *     E = 1/2 (lambda^2 - 1).
 *
 * The constitutive model returns the work-conjugate PK2 stress S. The truss
 * force relation
 *
 *     N = A0 lambda S
 *
 * is reported as the axial physical stress
 *
 *     sigma = lambda S.
 *
 * T3 currently supports reference recovery, exact recovery at the supplied
 * displacement and affine continuation about a supplied linearization state.
 */
void T3::compute_stress_strain(
    Field*           strain,
    Field*           stress,
    const Field&     target_displacement,
    const Field*     target_temperature,
    const RowMatrix& rst,
    const Field*     base_displacement,
    const Field*     base_temperature
) {
    // First element-nodal result row belonging to this element.
    Index offset = static_cast<Index>(this->elem_nodal_offset);

    // Validate all requirements before performing the actual recovery.
    logging::error(strain != nullptr || stress != nullptr,
        "T3: compute_stress_strain requires at least one output field");
    logging::error((!strain || strain->domain == FieldDomain::ELEMENT_NODAL)
                && (!stress || stress->domain == FieldDomain::ELEMENT_NODAL),
        "T3: stress/strain recovery requires ELEMENT_NODAL output");
    logging::error(rst.cols() >= 1,
        "T3: stress/strain coordinates require at least one natural coordinate");

    const auto [strain_value, stress_value] =
        evaluate_axial_response(
            target_displacement,
            target_temperature,
            base_displacement,
            base_temperature);

    // The axial strain and stress are constant over the two-node truss.
    for (Index i = 0; i < static_cast<Index>(rst.rows()); ++i) {
        const Index row = static_cast<Index>(offset) + i;

        if (strain) {
            for (Index component = 0; component < strain->components; ++component)
                (*strain)(row, component) = Precision(0);

            (*strain)(row, 0) = strain_value;
        }

        if (stress) {
            for (Index component = 0; component < stress->components; ++component)
                (*stress)(row, component) = Precision(0);

            (*stress)(row, 0) = stress_value;
        }
    }
}

/**
 * Recovers accumulated equivalent plastic strain from committed J2 history.
 *
 * T3 owns exactly one constitutive material point, hence PEEQ is constant over
 * the element and is copied to both element-nodal output rows. No constitutive
 * reevaluation is necessary: the accepted scalar history variable is read
 * directly from material_state_old.
 */
bool T3::compute_peeq(Field& peeq) {
    // First element-nodal result row belonging to this element.
    Index offset = static_cast<Index>(this->elem_nodal_offset);

    logging::error(peeq.domain == FieldDomain::ELEMENT_NODAL && peeq.components == 1,
        "T3: PEEQ recovery requires scalar ELEMENT_NODAL output");

    // Equivalent plastic strain is available only for the J2 material model.
    // Other elasticities simply do not contribute a PEEQ result.
    auto material = get_material();
    if (!material->has_elasticity()) return false;

    const auto* j2 = material->elasticity()->as<material::IsotropicJ2Elasticity>();
    if (!j2) return false;

    // T3 has one constitutive material point. Read the accumulated equivalent
    // plastic strain directly from the committed material history.
    Precision value = Precision(0);
    const auto& state = this->_model_data->material_state_old;

    if (state && state->components >= j2->state_size()) {
        const Precision* old_state = &(*state)(this->mp_index(0), 0);
        value = j2->equivalent_plastic_strain(old_state);
    }

    // The material-point value is constant over the truss and is therefore
    // written to both element-nodal result rows.
    peeq(static_cast<Index>(offset) + 0, 0) = value;
    peeq(static_cast<Index>(offset) + 1, 0) = value;

    return true;
}

/**
 * Computes the linearized element compliance contribution
 *
 *     J_e = u_e^T K_0 u_e.
 *
 * The stiffness request uses u0 = 0 and no force/geometric output, therefore
 * evaluate() returns the reference tangent only. This keeps compliance
 * consistent with the linear topology-optimization operator.
 */
void T3::compute_compliance(Field& displacement, Field& result) {
    // The compliance contribution of one linear truss element is
    //
    //     J_e = u_e^T K_0 u_e,
    //
    // where K_0 is the tangent evaluated at the reference state u0 = 0.
    Precision buffer[N * 3 * N * 3] {};
    const MapMatrix K = evaluate(
        buffer, nullptr, nullptr,
        nullptr, nullptr,
        nullptr, nullptr,
        false);

    // Gather the six translational element DOFs in local node-major ordering:
    //
    //     u_e = [u1_x, u1_y, u1_z, u2_x, u2_y, u2_z]^T.
    StaticVector<N * 3> u;
    for (Index node = 0; node < N; ++node) {
        u.template segment<3>(3 * node) = displacement.row_vec3(static_cast<Index>(node_ids[node]));
    }

    // Store the scalar element contribution
    //
    //     J_e = u_e^T K_0 u_e.
    result(static_cast<Index>(this->elem_id), 0) = u.dot(K * u);
}

/**
 * Recovers the axial section force from an exact base state followed by one
 * affine perturbation.
 *
 * The optional linearization field defines the exact base configuration u0.
 * Without one, u0 = 0 and the result is the ordinary reference-linearized
 * section force. For nonlinear output the caller supplies linearization =
 * displacement, so Delta u = 0 and the exact converged resultant is recovered.
 *
 * At the base state
 *
 *     lambda0 = l0 / L0,
 *     E0      = 1/2 (lambda0^2 - 1),
 *     S0      = S(E0 - E_th),
 *
 * while the displacement perturbation gives
 *
 *     Delta lambda = n0 . Delta r / L0,
 *     Delta E      = lambda0 Delta lambda,
 *     Delta S      = C0 Delta E.
 *
 * Linearizing
 *
 *     N = A0 lambda S
 *
 * about the base state yields
 *
 *     N ~= A0 [lambda0 S0 + S0 Delta lambda + lambda0 Delta S].
 *
 * This is identical to multiplying the physical axial stress recovered by
 * evaluate_axial_response() with the reference cross-sectional area A0.
 */
bool T3::compute_beam_section_forces(
    Field&       section_forces,
    const Field& target_displacement,
    const Field* target_temperature,
    const Field* base_displacement,
    const Field* base_temperature
) {
    // First element-nodal result row belonging to this element.
    Index offset = static_cast<Index>(this->elem_nodal_offset);

    const Precision stress_value =
        evaluate_axial_response(
            target_displacement,
            target_temperature,
            base_displacement,
            base_temperature).second;

    const Precision axial_force =
        get_section()->area_ * stress_value;

    // The axial section force is constant over the two-node truss.
    for (Index node = 0; node < N; ++node) {
        const Index row = static_cast<Index>(offset) + node;

        for (Index component = 0; component < section_forces.components; ++component) {
            section_forces(row, component) = Precision(0);
        }

        section_forces(row, 0) = axial_force;
    }

    return true;
}

} // namespace model
} // namespace fem