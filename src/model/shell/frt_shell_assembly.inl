/**
 * @file frt_shell_assembly.inl
 * @brief Implements shell-section evaluation and consistent element assembly.
 *
 * The generalized shell section maps local generalized strains to stress
 * resultants and a section tangent. The element routines integrate the
 * internal force, material tangent and geometric tangent over the actual
 * isoparametric reference midsurface.
 *
 * The physical geometric tangent is assembled by directly contracting the
 * generalized resultants with structured strain Hessian blocks. The objective
 * drilling stabilization is derived independently from one quadratic
 * in-plane-spin potential and contributes matching force, material-like and
 * geometric tangent terms.
 *
 * Generalized resultants and constitutive tangents remain local to the active
 * element evaluation. Only the physical nonlinear tangent path may write trial
 * material history; linear stiffness and prestress stiffness are state-neutral.
 *
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#include "frt_shell.h"

#include "../../core/logging.h"

#include <vector>

namespace fem::model {

/**
 * Evaluates shell-section resultants at every integration point.
 *
 * Generalized strains are supplied in the pointwise orthonormal reference
 * basis. The shell section receives the physical reference position and the
 * identical global basis used by the kinematic strain transformation. The
 * current consistent tangent is retained only when the caller requested B
 * matrices and therefore can consume a material tangent.
 *
 * Each in-plane integration point owns a contiguous block of section material
 * points. Every call reads the committed state. A writable trial-state pointer
 * is supplied only when `data.write_material_state` marks the physical nonlinear
 * equilibrium evaluation; state-neutral operators pass `nullptr` instead.
 *
 * @param data Active thread-local evaluation view.
 */
template<Index N>
void FRTShell<N>::compute_material_resultants(EvaluationData& data) const {
    logging::error(data.with_strain,
                   "FRTShell: material resultants require strain evaluation");

    ShellSection*   section      = shell_section();
    const Precision scale        = topology_stiffness_scale();
    const auto&     points       = reference_data().ip_points;
    const Index     state_stride = this->_model_data->material_state_old->components;

    // Evaluate one generalized section response for every shell integration point.
    for (Index ip = 0; ip < static_cast<Index>(points.size()); ++ip) {
        const std::size_t id = static_cast<std::size_t>(ip);
        const ReferencePoint& point = points[id];
        const Vec8& strain_values = data.ip_strain[id];

        ShellGeneralizedStrain strain(strain_values);
        ShellStressResultants  resultants;
        Mat8                   tangent;

        // Always read the committed through-thickness history. Only the physical
        // nonlinear tangent is allowed to write the corresponding trial rows.
        const Index      state_row = this->mp_index(ip, 0);
        const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
        Precision* new_state = data.write_material_state
            ? &(*this->_model_data->material_state_new)(state_row, 0)
            : nullptr;

        // Supply the same pointwise global basis used by the strain transformation.
        Mat3 basis = point.basis;

        section->evaluate(
            reference_position(point.r, point.s),
            basis,
            strain,
            old_state,
            new_state,
            state_stride,
            true,
            resultants,
            tangent
        );

        data.ip_resultants[id] = scale * resultants.values();

        // Retain the constitutive tangent only for paths that assemble B^T H B.
        if (!data.ip_tangent.empty()) {
            data.ip_tangent[id] = scale * tangent;
        }
    }
}

/**
 * Assembles the material part of the consistent element tangent.
 *
 * The integration follows
 *
 *     K_mat = integral_A0 B^T H B dA0,
 *
 * where `H` is the current generalized shell-section tangent and the physical
 * reference-area weight is evaluated directly from the cached reference point.
 *
 * @param data Active element evaluation data containing B and section tangents.
 * @param Kmat Output material tangent matrix.
 */
template<Index N>
void FRTShell<N>::assemble_material_stiffness(
    const EvaluationData& data,
    Mat6N&                Kmat
) const {
    logging::error(data.with_B,
                   "FRTShell: material stiffness requires B evaluation");

    Kmat.setZero();

    const auto& points = reference_data().ip_points;

    // Integrate the constitutive tangent over the actual curved reference area
    for (Index ip = 0; ip < static_cast<Index>(data.ip_B.size()); ++ip) {
        const std::size_t id = static_cast<std::size_t>(ip);
        const Precision weight = points[id].w * points[id].detJ;

        Kmat.noalias() += weight
                        * data.ip_B[id].transpose()
                        * data.ip_tangent[id]
                        * data.ip_B[id];
    }
}

/**
 * Adds one weighted compatible natural strain Hessian directly to the physical
 * geometric tangent.
 *
 * The eight generalized natural strain Hessians are never materialized. Metric
 * terms insert their constant translational identity blocks directly. Curvature
 * and shear terms use compact first and second nodal SO(3) derivatives for the
 * mixed and rotation-rotation blocks.
 *
 * @param data Active evaluation data containing second rotation derivatives.
 * @param point Integration or tying point whose compatible Hessian is weighted.
 * @param weights Generalized natural resultant weights.
 * @param Kgeo Element geometric tangent to update.
 */
template<Index N>
void FRTShell<N>::add_weighted_natural_hessian(
    const EvaluationData& data,
    const ReferencePoint& point,
    const Vec8&           weights,
    Mat6N&                Kgeo
) const {
    using Component = ShellGeneralizedStrain::Component;

    constexpr Index epsilon_rr = static_cast<Index>(Component::EpsilonXX);
    constexpr Index epsilon_ss = static_cast<Index>(Component::EpsilonYY);
    constexpr Index gamma_rs   = static_cast<Index>(Component::GammaXY);
    constexpr Index kappa_rr   = static_cast<Index>(Component::KappaXX);
    constexpr Index kappa_ss   = static_cast<Index>(Component::KappaYY);
    constexpr Index kappa_rs   = static_cast<Index>(Component::KappaXY);
    constexpr Index gamma_r3   = static_cast<Index>(Component::GammaXZ);
    constexpr Index gamma_s3   = static_cast<Index>(Component::GammaYZ);

    if (weights.cwiseAbs().maxCoeff() == Precision(0)) {
        return;
    }

    logging::error(data.rotations != nullptr,
        "FRTShell: geometric tangent requires nodal rotation derivatives");

    const auto& rotations = *data.rotations;
    const auto& ref       = reference_data();

    const VecN shape_r = point.shape_rs.col(0);
    const VecN shape_s = point.shape_rs.col(1);
    const VecN shape   = point.shape;

    // Membrane metric Hessians contain only direct translational identities
    add_xx_hessian<N>(shape_r, shape_r, Precision(0.5) * weights(epsilon_rr), Kgeo);
    add_xx_hessian<N>(shape_s, shape_s, Precision(0.5) * weights(epsilon_ss), Kgeo);
    add_xx_hessian<N>(shape_r, shape_s, weights(gamma_rs), Kgeo);

    // Curvature Hessians contain mixed translation-rotation and local
    // rotation-rotation SO(3) blocks
    add_xd_hessian<N>(data.state.x, rotations, ref.d0, shape_r, shape_r, weights(kappa_rr), Kgeo);
    add_xd_hessian<N>(data.state.x, rotations, ref.d0, shape_s, shape_s, weights(kappa_ss), Kgeo);
    add_xd_hessian<N>(data.state.x, rotations, ref.d0, shape_r, shape_s, weights(kappa_rs), Kgeo);
    add_xd_hessian<N>(data.state.x, rotations, ref.d0, shape_s, shape_r, weights(kappa_rs), Kgeo);

    // Transverse-shear Hessians use the interpolated current director field
    add_xd_hessian<N>(data.state.x, rotations, ref.d0, shape_r, shape, weights(gamma_r3), Kgeo);
    add_xd_hessian<N>(data.state.x, rotations, ref.d0, shape_s, shape, weights(gamma_s3), Kgeo);
}

/**
 * Assembles the physical stress-dependent geometric element tangent.
 *
 * Local resultants are multiplied by the complete reference-area weight, pulled
 * back through the pointwise local-basis transformation and then through the
 * transpose of the concrete MITC operator. The resulting compatible weights
 * act at the integration point itself and at all tying points.
 *
 * The tying buffer belongs to the active thread-local workspace and is reset in
 * place for every integration point, so geometric assembly performs no dynamic
 * allocation.
 *
 * @param data Active evaluation data containing resultants and second rotation
 * derivatives.
 * @param Kgeo Output physical geometric tangent matrix.
 */
template<Index N>
void FRTShell<N>::assemble_geometric_stiffness(
    const EvaluationData& data,
    Mat6N&                Kgeo
) const {
    logging::error(data.with_G,
        "FRTShell: geometric stiffness requires second rotation derivatives");
    logging::error(data.with_resultants,
        "FRTShell: geometric stiffness requires shell resultants");

    Kgeo.setZero();

    const auto& points = reference_data().ip_points;
    const auto& tying  = reference_data().tying_points;

    logging::error(data.geometric_tying_weights.size() == tying.size(),
                   "FRTShell: invalid geometric tying workspace size");

    for (Index ip = 0; ip < static_cast<Index>(points.size()); ++ip) {
        const std::size_t id = static_cast<std::size_t>(ip);
        const Precision weight = points[id].w * points[id].detJ;
        const Vec8 local_resultants = weight * data.ip_resultants[id];

        // Apply the transpose of the pointwise natural-to-local strain map
        // directly here because this pull-back is used only by geometric
        // assembly and does not justify a separate one-call helper.
        const Precision t00 = points[id].invJ(0, 0);
        const Precision t01 = points[id].invJ(0, 1);
        const Precision t10 = points[id].invJ(1, 0);
        const Precision t11 = points[id].invJ(1, 1);

        StaticMatrix<3, 3> in_plane;
        in_plane << t00 * t00,               t01 * t01,               t00 * t01,
                    t10 * t10,               t11 * t11,               t10 * t11,
                    Precision(2) * t00 * t10, Precision(2) * t01 * t11,
                    t00 * t11 + t01 * t10;

        Vec8 natural_resultants = Vec8::Zero();
        natural_resultants.template segment<3>(0) = in_plane.transpose() * local_resultants.template segment<3>(0);
        natural_resultants.template segment<3>(3) = in_plane.transpose() * local_resultants.template segment<3>(3);
        natural_resultants.template segment<2>(6) = points[id].invJ.transpose() * local_resultants.template segment<2>(6);

        Vec8 compatible_weights = Vec8::Zero();
        for (Vec8& tying_weight : data.geometric_tying_weights) {
            tying_weight.setZero();
        }

        // Apply the exact transpose of the topology-specific MITC interpolation
        pull_back_mitc_resultants(
            points[id],
            natural_resultants,
            compatible_weights,
            data.geometric_tying_weights
        );

        // Add the compatible integration-point Hessian contribution
        add_weighted_natural_hessian(data, points[id], compatible_weights, Kgeo);

        // Add all compatible tying-point Hessian contributions
        for (Index tying_id = 0; tying_id < static_cast<Index>(tying.size()); ++tying_id) {
            add_weighted_natural_hessian(
                data,
                tying[static_cast<std::size_t>(tying_id)],
                data.geometric_tying_weights[static_cast<std::size_t>(tying_id)],
                Kgeo
            );
        }
    }

    // Remove only round-off asymmetry from the analytically symmetric tangent
    Kgeo = Precision(0.5) * (Kgeo + Kgeo.transpose());
}

/**
 * Assembles the nonlinear physical shell internal force vector.
 *
 * The generalized shell resultants are integrated through
 *
 *     f_int = integral_A0 B^T n dA0.
 *
 * @param data Active evaluation data containing B matrices and resultants.
 * @param internal_force Output element internal force vector.
 */
template<Index N>
void FRTShell<N>::assemble_internal_force(
    const EvaluationData& data,
    Vec6N&                internal_force
) const {
    logging::error(data.with_B,
                   "FRTShell: internal force requires B evaluation");
    logging::error(data.with_resultants,
                   "FRTShell: internal force requires shell resultants");

    internal_force.setZero();

    const auto& points = reference_data().ip_points;

    // Integrate the physical shell resultants over the curved reference area
    for (Index ip = 0; ip < static_cast<Index>(data.ip_B.size()); ++ip) {
        const std::size_t id = static_cast<std::size_t>(ip);
        const Precision weight = points[id].w * points[id].detJ;

        internal_force.noalias() +=
            weight * data.ip_B[id].transpose() * data.ip_resultants[id];
    }
}

/**
 * Assembles objective drilling stabilization from one quadratic potential.
 *
 * At every integration point the interpolated nodal rotation field acts on the
 * two pointwise reference tangents,
 *
 *     a1 = sum_i N_i R_i e1,
 *     a2 = sum_i N_i R_i e2,
 *
 * while `x_,a` and `x_,b` are the current midsurface tangents measured with
 * respect to the same orthonormal reference coordinates. The drilling strain is
 *
 *     gamma_d = 1/2 (a1 . x_,b - a2 . x_,a).
 *
 * Under an arbitrary finite rigid-body rotation `Q`, all four vectors rotate by
 * `Q`, so `gamma_d` remains exactly zero. Its small-strain limit is
 *
 *     gamma_d = theta_3 - 1/2 (u_2,1 - u_1,2),
 *
 * which is the standard difference between the independent drilling rotation
 * and the in-plane continuum spin.
 *
 * The stabilization potential is
 *
 *     Pi_d = 1/2 integral_A0 k_d gamma_d^2 dA0,
 *
 * with `k_d = drill_scale * |A66|` evaluated from the zero-strain shell-section
 * tangent. The force and tangent are therefore
 *
 *     f_d = integral_A0 k_d gamma_d B_d^T dA0,
 *     K_d = integral_A0 k_d (B_d^T B_d + gamma_d G_d) dA0.
 *
 * Both quantities are assembled directly. No drilling B or Hessian fields are
 * retained in the thread-local workspace.
 *
 * @param data Active evaluation data containing compact SO(3) derivatives.
 * @param stiffness_matrix Optional tangent matrix to update.
 * @param internal_force Optional internal force vector to update.
 */
template<Index N>
void FRTShell<N>::assemble_drill_stabilization(
    const EvaluationData& data,
    Mat6N*                stiffness_matrix,
    Vec6N*                internal_force
) const {
    logging::error(data.with_B,
        "FRTShell: drilling stabilization requires first rotation derivatives");
    logging::error(data.rotations != nullptr,
        "FRTShell: drilling stabilization requires nodal rotations");
    logging::error(data.ip_drill_stiffness.size() == reference_data().ip_points.size(),
        "FRTShell: invalid drilling stiffness workspace size");

    const auto& rotations = *data.rotations;
    const auto& points    = reference_data().ip_points;

    for (Index ip = 0; ip < static_cast<Index>(points.size()); ++ip) {
        const std::size_t id = static_cast<std::size_t>(ip);
        const ReferencePoint& point = points[id];
        const Precision k_d = data.ip_drill_stiffness[id];

        if (k_d == Precision(0)) {
            continue;
        }

        Vec3 x_a = Vec3::Zero();
        Vec3 x_b = Vec3::Zero();
        Vec3 a1  = Vec3::Zero();
        Vec3 a2  = Vec3::Zero();

        // Interpolate current surface tangents and the two independently
        // rotated reference tangent vectors
        for (Index node = 0; node < num_nodes; ++node) {
            const Vec3 x_i = data.state.x.row(node).transpose();
            x_a += point.shape_ab.col(0)(node) * x_i;
            x_b += point.shape_ab.col(1)(node) * x_i;
            a1  += point.shape(node) * rotations[node].value * point.basis.col(0);
            a2  += point.shape(node) * rotations[node].value * point.basis.col(1);
        }

        const Precision gamma_d = Precision(0.5)
                                * (a1.dot(x_b) - a2.dot(x_a));

        Vec6N B_d = Vec6N::Zero();

        // Translational derivative:
        // d(gamma_d)/du_i = 1/2 (N_i,b a1 - N_i,a a2)
        for (Index node = 0; node < num_nodes; ++node) {
            const Index base = dofs_per_node * node;
            const Vec3 derivative = Precision(0.5)
                                  * (point.shape_ab.col(1)(node) * a1
                                   - point.shape_ab.col(0)(node) * a2);
            B_d.template segment<3>(base) = derivative;
        }

        // Rotational derivative:
        // d(gamma_d)/dtheta_ia = N_i/2 [dR_ia e1 . x_b - dR_ia e2 . x_a]
        for (Index node = 0; node < num_nodes; ++node) {
            const Index rot_base = dofs_per_node * node + 3;
            const Precision shape = point.shape(node);

            for (Index a = 0; a < 3; ++a) {
                const Vec3 da1 = shape * rotations[node].d1[a] * point.basis.col(0);
                const Vec3 da2 = shape * rotations[node].d1[a] * point.basis.col(1);
                B_d(rot_base + a) = Precision(0.5)
                                  * (da1.dot(x_b) - da2.dot(x_a));
            }
        }

        const Precision weighted_stiffness = point.w * point.detJ * k_d;

        // Add the first variation of the quadratic drilling potential
        if (internal_force) {
            internal_force->noalias() += weighted_stiffness * gamma_d * B_d;
        }

        if (!stiffness_matrix) {
            continue;
        }

        // Add the positive semidefinite B_d^T B_d contribution
        stiffness_matrix->noalias() += weighted_stiffness * B_d * B_d.transpose();

        // The second kinematic derivative is required only for the complete
        // nonlinear tangent. Linear shell stiffness at the reference state has
        // gamma_d = 0, so this term vanishes there identically.
        if (!data.with_G || gamma_d == Precision(0)) {
            continue;
        }

        const Precision geometric_scale = weighted_stiffness * gamma_d;

        // Mixed translation-rotation blocks arise because the rotated tangent
        // vectors multiply the current midsurface derivatives
        for (Index rot_node = 0; rot_node < num_nodes; ++rot_node) {
            const Index rot_base = dofs_per_node * rot_node + 3;
            const Precision shape = point.shape(rot_node);

            for (Index a = 0; a < 3; ++a) {
                const Vec3 da1 = shape * rotations[rot_node].d1[a] * point.basis.col(0);
                const Vec3 da2 = shape * rotations[rot_node].d1[a] * point.basis.col(1);

                for (Index x_node = 0; x_node < num_nodes; ++x_node) {
                    const Index x_base = dofs_per_node * x_node;
                    const Vec3 mixed = Precision(0.5)
                                     * (point.shape_ab.col(1)(x_node) * da1
                                      - point.shape_ab.col(0)(x_node) * da2);

                    stiffness_matrix->template block<3, 1>(x_base, rot_base + a)
                        += geometric_scale * mixed;
                    stiffness_matrix->template block<1, 3>(rot_base + a, x_base)
                        += geometric_scale * mixed.transpose();
                }
            }

            // Pure rotation-rotation terms remain local to one nodal SO(3)
            // parameterization because the interpolated rotation field is a
            // linear sum of independent nodal rotation matrices
            for (Index a = 0; a < 3; ++a) {
                for (Index b = 0; b < 3; ++b) {
                    const Vec3 d2a1 = shape * rotations[rot_node].d2[a][b] * point.basis.col(0);
                    const Vec3 d2a2 = shape * rotations[rot_node].d2[a][b] * point.basis.col(1);
                    const Precision second = Precision(0.5)
                                           * (d2a1.dot(x_b) - d2a2.dot(x_a));

                    (*stiffness_matrix)(rot_base + a, rot_base + b) +=
                        geometric_scale * second;
                }
            }
        }
    }
}

/**
 * Assembles the linear/reference shell stiffness.
 *
 * The undeformed shell state supplies the reference B matrices. Section
 * tangents are queried at zero generalized strain without a writable material
 * target state, and the reference drilling tangent is added from the same
 * objective stabilization potential. No current resultants or geometric
 * stiffness enter this operator.
 *
 * @param buffer Caller-provided dense element matrix storage.
 * @return Mapped linear shell stiffness.
 */
template<Index N>
MapMatrix FRTShell<N>::stiffness(Precision* buffer) {
    const CurrentState state = reference_state();
    const EvaluationData data = init_evaluation(
        state,
        true,
        true,
        false,
        false,
        false
    );

    Mat6N Kmat;
    assemble_material_stiffness(data, Kmat);
    assemble_drill_stabilization(data, &Kmat, nullptr);

    // Remove only numerical asymmetry from the analytically symmetric tangent.
    Kmat = Precision(0.5) * (Kmat + Kmat.transpose());

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = Kmat;
    return mapped;
}

/**
 * Evaluates the physical geometric stiffness from a linearized prestress state.
 *
 * Linear buckling supplies the displacement obtained from the linear preload
 * solve. Prestress must therefore use the same reference B matrices and the
 * linearized shell constitutive measure rather than finite-rotation shell
 * strains. At every integration point the element evaluates
 *
 *     epsilon_0 = B_0 u,
 *     n_0       = n(epsilon_0, state_committed),
 *
 * with `use_green_lagrange = false`. The resulting generalized resultants are
 * contracted with the strain Hessians evaluated at the reference state. The
 * complete path is state-neutral and keeps resultants local to this element.
 *
 * @param buffer Caller-provided dense element matrix storage.
 * @param displacement Global nodal displacement field from the linear preload solve.
 * @return Mapped physical geometric stiffness matrix.
 */
template<Index N>
MapMatrix FRTShell<N>::stiffness_geom(
    Precision*   buffer,
    const Field& displacement
) {
    const CurrentState state = reference_state();
    EvaluationData data = init_evaluation(
        state,
        true,
        true,
        true,
        false,
        false
    );

    const Vec6N q = element_displacement_vector(displacement);
    const auto& points = reference_data().ip_points;
    const Index state_stride = this->_model_data->material_state_old->components;
    const Precision scale = topology_stiffness_scale();
    ShellSection* section = shell_section();

    // Keep the prestress resultants local and reuse the allocation across calls
    // on the same worker thread. They are attached to EvaluationData only for
    // the subsequent geometric contraction.
    thread_local std::vector<Vec8> linear_resultants;
    linear_resultants.resize(points.size());

    for (Index ip = 0; ip < static_cast<Index>(points.size()); ++ip) {
        const std::size_t id = static_cast<std::size_t>(ip);
        const ReferencePoint& point = points[id];
        const ShellGeneralizedStrain strain(data.ip_B[id] * q);
        ShellStressResultants resultants;
        Mat8 tangent;

        const Index      state_row = this->mp_index(ip, 0);
        const Precision* old_state = &(*this->_model_data->material_state_old)(state_row, 0);
        Mat3 basis = point.basis;

        section->evaluate(
            reference_position(point.r, point.s),
            basis,
            strain,
            old_state,
            nullptr,
            state_stride,
            false,
            resultants,
            tangent
        );

        linear_resultants[id] = scale * resultants.values();
    }

    data.ip_resultants      = Span<Vec8>(linear_resultants);
    data.with_resultants    = true;

    Mat6N Kgeo;
    assemble_geometric_stiffness(data, Kgeo);

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = Kgeo;
    return mapped;
}

/**
 * Assembles nonlinear shell internal force and optionally the consistent tangent.
 *
 * The supplied displacement field defines the physical trial configuration.
 * Every integration point evaluates its generalized strain and section response
 * once from committed material history into the persistent trial state. The
 * resulting local resultants are used immediately for the internal force and,
 * when matrix storage is requested, for the material and geometric tangent.
 *
 * A null matrix buffer performs a residual-only evaluation. It still updates
 * the physical material trial state and assembles the complete internal force,
 * including objective drilling stabilization, but skips all tangent-only work.
 *
 * @param buffer Optional dense tangent storage; null requests force only.
 * @param nodal_forces Global nodal internal-force field to increment.
 * @param displacement Trial nodal displacement field.
 * @return Mapped complete tangent, or an empty map for residual-only evaluation.
 */
template<Index N>
MapMatrix FRTShell<N>::stiffness_tangent(
    Precision*   buffer,
    NodeData&    nodal_forces,
    const Field& displacement
) {
    logging::error(nodal_forces.components >= dofs_per_node,
                   "FRTShell: nonlinear internal force requires six nodal components");

    const bool with_tangent = buffer != nullptr;
    const CurrentState state = current_state_from_displacement(displacement);
    const EvaluationData data = init_evaluation(
        state,
        true,
        true,
        with_tangent,
        true,
        true
    );

    Vec6N internal_force;
    assemble_internal_force(data, internal_force);

    Mat6N tangent;
    if (with_tangent) {
        Mat6N Kmat;
        Mat6N Kgeo;
        assemble_material_stiffness(data, Kmat);
        assemble_geometric_stiffness(data, Kgeo);
        tangent = Kmat + Kgeo;
        assemble_drill_stabilization(data, &tangent, &internal_force);
        tangent = Precision(0.5) * (tangent + tangent.transpose());
    } else {
        assemble_drill_stabilization(data, nullptr, &internal_force);
    }

    // Scatter the complete element internal force into the global nodal field.
    for (Index node = 0; node < num_nodes; ++node) {
        const Index node_id = static_cast<Index>(this->node_ids[node]);
        for (Index dof = 0; dof < dofs_per_node; ++dof) {
            nodal_forces(node_id, dof) += internal_force(dofs_per_node * node + dof);
        }
    }

    if (!with_tangent) {
        return MapMatrix(nullptr, 0, 0);
    }

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = tangent;
    return mapped;
}

} // namespace fem::model
