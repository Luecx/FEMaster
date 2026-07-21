/**
 * @file frt_shell_assembly.cpp
 * @brief Implements shell-section evaluation and consistent element assembly.
 *
 * The generalized shell section maps local generalized strains to stress
 * resultants and a section tangent. The element routines integrate the
 * internal force, material tangent and geometric tangent over the actual
 * isoparametric reference midsurface.
 *
 * A small consistent drilling contribution regularizes rotation about the
 * nodal reference director. The drilling force and stiffness are derived from
 * the same quadratic penalty potential.
 *
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 20.07.2026
 */

#include "frt_shell.h"

#include "../../core/logging.h"

#include <cmath>

namespace fem::model {

/**
 * Loads previously stored generalized shell resultants at the integration
 * points.
 *
 * The field components follow the ordering
 * `[N11,N22,N12,M11,M22,M12,Q13,Q23]`.
 */
template<Index N>
void FRTShell<N>::load_ip_resultants(EvaluationData& data,
                                     const Field&    ip_stress,
                                     int             ip_start_idx) const {
    logging::error(ip_stress.components >= num_strains,
                   "FRTShell: shell resultant field requires eight components "
                   "[N11,N22,N12,M11,M22,M12,Q13,Q23]");

    for (Index ip = 0; ip < static_cast<Index>(data.ip_resultants.size()); ++ip) {
        const Index row = static_cast<Index>(ip_start_idx) + ip;

        for (Index component = 0; component < num_strains; ++component) {
            data.ip_resultants[static_cast<std::size_t>(ip)](component) = ip_stress(row, component);
        }
    }
}

/**
 * Evaluates the shell section at every integration point.
 *
 * The generalized strains are expressed in the pointwise orthonormal reference
 * basis. The section receives the physical reference position and the same
 * local basis used by the kinematic strain vector.
 */
template<Index N>
void FRTShell<N>::compute_material_resultants(EvaluationData& data) const {
    logging::error(data.with_strain,
                   "FRTShell: material resultants require strain evaluation");

    ShellSection*   section = shell_section();
    const Precision scale   = topology_stiffness_scale();
    const auto&     points  = reference_data().ip_points;

    for (Index ip = 0; ip < static_cast<Index>(points.size()); ++ip) {
        const ReferencePoint& point = points[static_cast<std::size_t>(ip)];
        const Vec8& strain_values = data.ip_strain[static_cast<std::size_t>(ip)];

        ShellGeneralizedStrain strain(strain_values);
        ShellStressResultants  resultants;
        Mat8                   tangent;

        Mat3 basis;
        basis.col(0) = point.e1;
        basis.col(1) = point.e2;
        basis.col(2) = point.e3;

        section->evaluate(
            reference_position(point.r, point.s),
            basis,
            strain,
            true,
            resultants,
            tangent
        );

        data.ip_resultants[static_cast<std::size_t>(ip)] = scale * resultants.values();
        data.ip_tangent[static_cast<std::size_t>(ip)]    = scale * tangent;
    }
}

/**
 * Assembles the material part of the consistent element tangent.
 *
 * The integration follows
 *
 *     K_mat = integral_A0 B^T H B dA0,
 *
 * where `H` is the generalized shell-section tangent and every integration
 * weight already contains the complete reference surface Jacobian.
 */
template<Index N>
void FRTShell<N>::assemble_material_stiffness(const EvaluationData& data,
                                              Mat6N&                Kmat) const {
    logging::error(data.with_B,
                   "FRTShell: material stiffness requires B evaluation");

    Kmat.setZero();

    for (Index ip = 0; ip < static_cast<Index>(data.ip_B.size()); ++ip) {
        const std::size_t id = static_cast<std::size_t>(ip);
        Kmat.noalias() += data.ip_weight[id]
                        * data.ip_B[id].transpose()
                        * data.ip_tangent[id]
                        * data.ip_B[id];
    }
}

/**
 * Assembles the geometric part of the consistent element tangent.
 *
 * Every generalized strain component contributes its exact Hessian weighted by
 * the associated generalized stress resultant:
 *
 *     K_geo = integral_A0 sum_a n[a] G[a] dA0.
 */
template<Index N>
void FRTShell<N>::assemble_geometric_stiffness(const EvaluationData& data,
                                               Mat6N&                Kgeo) const {
    logging::error(data.with_G,
                   "FRTShell: geometric stiffness requires G evaluation");
    logging::error(data.with_resultants,
                   "FRTShell: geometric stiffness requires shell resultants");

    Kgeo.setZero();

    for (Index ip = 0; ip < static_cast<Index>(data.ip_G.size()); ++ip) {
        const std::size_t id = static_cast<std::size_t>(ip);

        for (Index component = 0; component < num_strains; ++component) {
            Kgeo.noalias() += data.ip_weight[id]
                            * data.ip_resultants[id](component)
                            * data.ip_G[id][component];
        }
    }

    // Remove round-off asymmetry from the analytically symmetric tangent
    Kgeo = Precision(0.5) * (Kgeo + Kgeo.transpose());
}

/**
 * Assembles the nonlinear internal force vector.
 *
 * The generalized stress resultants are integrated through
 *
 *     f_int = integral_A0 B^T n dA0.
 */
template<Index N>
void FRTShell<N>::assemble_internal_force(const EvaluationData& data,
                                          Vec6N&                internal_force) const {
    logging::error(data.with_B,
                   "FRTShell: internal force requires B evaluation");
    logging::error(data.with_resultants,
                   "FRTShell: internal force requires shell resultants");

    internal_force.setZero();

    for (Index ip = 0; ip < static_cast<Index>(data.ip_B.size()); ++ip) {
        const std::size_t id = static_cast<std::size_t>(ip);
        internal_force.noalias() += data.ip_weight[id]
                                  * data.ip_B[id].transpose()
                                  * data.ip_resultants[id];
    }
}

template<Index N>
Precision FRTShell<N>::drill_stiffness_per_node(const Mat8& H) const {
    const Precision shear_scale = Precision(0.5)
                                * (std::abs(H(6, 6)) + std::abs(H(7, 7)));

    return drill_scale
         * shear_scale
         * reference_data().area
         / Precision(num_nodes);
}

/**
 * Adds the tangent of the quadratic drilling penalty.
 *
 * The regularized rotational component is the projection of the total nodal
 * rotation vector onto the nodal reference director.
 */
template<Index N>
void FRTShell<N>::add_drill_stiffness(Mat6N& stiffness_matrix,
                                      const Mat8& H) const {
    const auto&     ref = reference_data();
    const Precision k   = drill_stiffness_per_node(H);

    for (Index node = 0; node < num_nodes; ++node) {
        const Vec3 d0 = ref.d0.row(node).transpose();

        stiffness_matrix.template block<3, 3>(
            dofs_per_node * node + 3,
            dofs_per_node * node + 3
        ) += k * d0 * d0.transpose();
    }
}

/**
 * Adds the internal force generated by the quadratic drilling penalty.
 */
template<Index N>
void FRTShell<N>::add_drill_force(Vec6N&             internal_force,
                                  const CurrentState& state,
                                  const Mat8&         H) const {
    const auto&     ref = reference_data();
    const Precision k   = drill_stiffness_per_node(H);

    for (Index node = 0; node < num_nodes; ++node) {
        const Vec3 d0    = ref.d0.row(node).transpose();
        const Vec3 theta = state.theta.row(node).transpose();
        const Precision drill_rotation = theta.dot(d0);

        internal_force.template segment<3>(dofs_per_node * node + 3) +=
            k * drill_rotation * d0;
    }
}

/**
 * Returns the material stiffness in the current element configuration.
 *
 * This interface is used by linear analyses and by algorithms that request the
 * constitutive part separately. The compatible/MITC B matrix is evaluated at
 * the current nodal state, while geometric stress-resultant terms are omitted.
 */
template<Index N>
MapMatrix FRTShell<N>::stiffness(Precision* buffer) {
    const CurrentState state = current_state();
    const EvaluationData data = init_evaluation(
        state,
        true,
        true,
        false,
        true
    );

    Mat6N Kmat;
    assemble_material_stiffness(data, Kmat);

    add_drill_stiffness(Kmat, data.H);
    Kmat = Precision(0.5) * (Kmat + Kmat.transpose());

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = Kmat;
    return mapped;
}

/**
 * Returns the geometric stiffness generated by a supplied integration-point
 * shell-resultant field.
 */
template<Index N>
MapMatrix FRTShell<N>::stiffness_geom(Precision*   buffer,
                                      const Field& ip_stress,
                                      int          ip_start_idx) {
    const CurrentState state = current_state();
    const EvaluationData data = init_evaluation(
        state,
        true,
        true,
        true,
        true,
        &ip_stress,
        ip_start_idx
    );

    Mat6N Kgeo;
    assemble_geometric_stiffness(data, Kgeo);

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = Kgeo;
    return mapped;
}

/**
 * Assembles the complete consistent nonlinear tangent and internal force.
 *
 * The supplied displacement field defines the trial configuration. Material
 * resultants are evaluated from the generalized strains, and both material and
 * geometric tangent contributions are assembled. The nodal internal-force
 * field is incremented in global degrees-of-freedom ordering.
 */
template<Index N>
MapMatrix FRTShell<N>::stiffness_tangent(Precision*   buffer,
                                         Field&       ip_stress_state,
                                         NodeData&    nodal_forces,
                                         const Field& displacement) {
    (void) ip_stress_state;

    const CurrentState state = current_state_from_displacement(displacement);
    const EvaluationData data = init_evaluation(
        state,
        true,
        true,
        true,
        true
    );

    Mat6N Kmat;
    Mat6N Kgeo;
    Vec6N internal_force;

    assemble_material_stiffness (data, Kmat);
    assemble_geometric_stiffness(data, Kgeo);
    assemble_internal_force     (data, internal_force);

    Mat6N tangent = Kmat + Kgeo;

    add_drill_stiffness(tangent, data.H);
    add_drill_force(internal_force, state, data.H);

    tangent = Precision(0.5) * (tangent + tangent.transpose());

    // Scatter the element internal force into the global nodal field
    for (Index node = 0; node < num_nodes; ++node) {
        const Index node_id = static_cast<Index>(this->node_ids[node]);

        for (Index dof = 0; dof < dofs_per_node; ++dof) {
            nodal_forces(node_id, dof) += internal_force(dofs_per_node * node + dof);
        }
    }

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = tangent;
    return mapped;
}

/**
 * Recovers the internal force from a stored shell-resultant field.
 *
 * This path is used when the material update has already been performed and
 * only `B^T n` must be integrated in the current configuration.
 */
template<Index N>
void FRTShell<N>::compute_internal_force_nonlinear(Field&       node_forces,
                                                   const Field& ip_stress) {
    const CurrentState state = current_state();
    const EvaluationData data = init_evaluation(
        state,
        true,
        true,
        false,
        true,
        &ip_stress,
        static_cast<int>(this->ip_index(0))
    );

    Vec6N internal_force;
    assemble_internal_force(data, internal_force);
    add_drill_force(internal_force, state, data.H);

    for (Index node = 0; node < num_nodes; ++node) {
        const Index node_id = static_cast<Index>(this->node_ids[node]);

        for (Index dof = 0; dof < dofs_per_node; ++dof) {
            node_forces(node_id, dof) += internal_force(dofs_per_node * node + dof);
        }
    }
}

} // namespace fem::model
