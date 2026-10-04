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
 * T3 currently supports only reference recovery or exact recovery at the
 * supplied displacement; intermediate affine recovery points are rejected.
 */
void T3::compute_stress_strain(
    Field*           strain,
    Field*           stress,
    const Field&     displacement,
    const RowMatrix& rst,
    int              offset,
    const Field*     linearization,
    const Field*     thermal_free_strain
) {
    (void) thermal_free_strain;

    logging::error(strain != nullptr || stress != nullptr,
        "T3: compute_stress_strain requires at least one output field");
    logging::error(rst.cols() >= 1,
        "T3: stress/strain coordinates require at least one natural coordinate");
    logging::error(linearization == nullptr || linearization == &displacement,
        "T3: intermediate recovery expansion points are not supported");

    const Precision L0 = length_reference();
    logging::error(L0 > Precision(0),
        "T3: zero reference length in compute_stress_strain for element ", this->elem_id);

    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
    auto             elasticity = get_elasticity();

    Precision strain_value = Precision(0);
    Precision stress_value = Precision(0);

    if (linearization == nullptr) {
        // ---------------------------------------------------------------------
        // infinitesimal reference recovery
        // ---------------------------------------------------------------------

        logging::error(elasticity->supports_axial_linearized(),
            "T3: material does not support linearized axial evaluation for element ",
            this->elem_id);

        const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
        const Vec3 u2 = displacement.row_vec3(static_cast<Index>(node_ids[1]));
        const Vec3 N0 = direction_reference();

        // Project the relative displacement onto the reference axis:
        //
        //     epsilon = N0 . (u2 - u1) / L0.
        const AxialStrainLinearized axial_strain(N0.dot(u2 - u1) / L0);
        AxialStressCauchy           axial_stress;

        elasticity->evaluate(
            axial_strain,
            old_state,
            nullptr,
            axial_stress,
            nullptr
        );

        strain_value = axial_strain.value();
        stress_value = axial_stress.value();

    } else {
        // ---------------------------------------------------------------------
        // exact Total-Lagrangian recovery
        // ---------------------------------------------------------------------

        logging::error(elasticity->supports_axial_green_lagrange(),
            "T3: material does not support Green-Lagrange axial evaluation for element ",
            this->elem_id);

        const Vec3 X1 = node_position_reference(0);
        const Vec3 X2 = node_position_reference(1);
        const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
        const Vec3 u2 = displacement.row_vec3(static_cast<Index>(node_ids[1]));

        // Build the exact requested configuration rather than relying on the
        // model's current POSITION field:
        //
        //     r = (X2 + u2) - (X1 + u1).
        const Vec3      axis   = (X2 + u2) - (X1 + u1);
        const Precision length = axis.norm();

        logging::error(length > Precision(0),
            "T3: zero length in exact stress recovery for element ", this->elem_id);

        const Precision lambda = length / L0;

        // The Green-Lagrange strain is constant over a two-node truss:
        //
        //     E = 1/2 (lambda^2 - 1).
        const AxialStrainGreenLagrange axial_strain =
            AxialStrainGreenLagrange::from_stretch(lambda);
        AxialStressPK2 axial_stress;

        elasticity->evaluate(
            axial_strain,
            old_state,
            nullptr,
            axial_stress,
            nullptr
        );

        strain_value = axial_strain.value();

        // For the adopted truss force measure,
        //
        //     N = A0 lambda S = A0 sigma,
        //
        // so the physical axial stress written to output is sigma = lambda S.
        stress_value = lambda * axial_stress.value();
    }

    // The axial strain/stress state is constant over the element. Replicate it
    // to every requested recovery coordinate and clear unsupported components.
    for (Index i = 0; i < static_cast<Index>(rst.rows()); ++i) {
        const Index row = static_cast<Index>(offset) + i;

        if (strain) {
            for (Index component = 0; component < strain->components; ++component) {
                (*strain)(row, component) = Precision(0);
            }
            (*strain)(row, 0) = strain_value;
        }

        if (stress) {
            for (Index component = 0; component < stress->components; ++component) {
                (*stress)(row, component) = Precision(0);
            }
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
bool T3::compute_peeq(Field& peeq, int offset) {
    logging::error(peeq.domain == FieldDomain::ELEMENT_NODAL && peeq.components == 1,
        "T3: PEEQ recovery requires scalar ELEMENT_NODAL output");

    auto material = get_material();
    if (!material || !material->has_elasticity()) {
        return false;
    }

    const auto* j2 = material->elasticity()->as<material::IsotropicJ2Elasticity>();
    if (!j2) {
        return false;
    }

    Precision value = Precision(0);

    const auto& state = this->_model_data->material_state_old;
    if (state && state->components >= j2->state_size()) {
        const Precision* old_state = &(*state)(this->mp_index(0), 0);
        value = j2->equivalent_plastic_strain(old_state);
    }

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
    Precision buffer[N * 3 * N * 3] {};
    MapMatrix K = evaluate(buffer, nullptr, nullptr, nullptr, nullptr, nullptr, false);

    StaticVector<N * 3> u;
    for (Index node = 0; node < N; ++node) {
        u.template segment<3>(3 * node) =
            displacement.row_vec3(static_cast<Index>(node_ids[node]));
    }

    result(static_cast<Index>(this->elem_id), 0) = u.dot(K * u);
}

/**
 * Recovers the linearized axial section force for beam-style output.
 *
 * The infinitesimal axial strain is
 *
 *     epsilon = N0 . (u2 - u1) / L0,
 *
 * the material returns Cauchy stress sigma, and the constant axial resultant is
 *
 *     N = A0 sigma.
 *
 * Only stress is required, so the constitutive tangent is intentionally omitted.
 */
bool T3::compute_beam_section_forces(
    Field&       section_forces,
    const Field& displacement,
    int          offset
) {
    const Precision L0 = length_reference();
    logging::error(L0 > Precision(0),
        "T3: zero reference length in compute_beam_section_forces for element ",
        this->elem_id);

    const Vec3 u1 = displacement.row_vec3(static_cast<Index>(node_ids[0]));
    const Vec3 u2 = displacement.row_vec3(static_cast<Index>(node_ids[1]));
    const Vec3 N0 = direction_reference();

    const AxialStrainLinearized axial_strain(N0.dot(u2 - u1) / L0);
    AxialStressCauchy           axial_stress;

    auto elasticity = get_elasticity();
    logging::error(elasticity->supports_axial_linearized(),
        "T3: material does not support linearized axial evaluation for element ",
        this->elem_id);

    const Index      state_row = this->mp_index(0);
    const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);

    elasticity->evaluate(
        axial_strain,
        old_state,
        nullptr,
        axial_stress,
        nullptr
    );

    const Precision axial_force = get_section()->area_ * axial_stress.value();

    // T3 carries only one constant axial resultant. All bending, shear and
    // torsional resultants are zero by construction.
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
