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
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#include "frt_shell.h"

#include "../../core/logging.h"

namespace fem::model {

namespace {

/**
 * Adds one stress-weighted Hessian of a scalar product between two current
 * midsurface derivatives directly to the translational element block.
 *
 * For
 *
 *     x_a = sum_i a_i x_i,
 *     x_b = sum_i b_i x_i,
 *
 * the second derivative of `scale * x_a . x_b` consists only of constant
 * translational identity blocks
 *
 *     scale * (a_i b_j + b_i a_j) I.
 *
 * The identities are inserted through the block diagonals and are never stored
 * as explicit derivative fields.
 *
 * @tparam N Number of shell nodes.
 * @param a_coefficients Interpolation coefficients of the first position field.
 * @param b_coefficients Interpolation coefficients of the second position field.
 * @param scale Complete stress and quadrature multiplier.
 * @param Kgeo Element geometric tangent to update.
 */
template<Index N>
void add_xx_weighted_hessian(
    const typename FRTShell<N>::VecN& a_coefficients,
    const typename FRTShell<N>::VecN& b_coefficients,
    Precision                         scale,
    typename FRTShell<N>::Mat6N&      Kgeo
) {
    if (scale == Precision(0)) {
        return;
    }

    using Shell = FRTShell<N>;

    // Insert every translational identity block directly into the global
    // element matrix without constructing an intermediate Hessian
    for (Index i = 0; i < N; ++i) {
        const Index row = Shell::dofs_per_node * i;

        for (Index j = 0; j < N; ++j) {
            const Index col = Shell::dofs_per_node * j;
            const Precision block_scale = scale
                * (a_coefficients(i) * b_coefficients(j)
                 + b_coefficients(i) * a_coefficients(j));

            if (block_scale != Precision(0)) {
                Kgeo.template block<3, 3>(row, col).diagonal().array() += block_scale;
            }
        }
    }
}

/**
 * Adds one stress-weighted Hessian of a scalar product between a current
 * midsurface derivative and a rotated director field.
 *
 * The scalar kinematic quantity is
 *
 *     scale * (sum_i a_i x_i) . (sum_j b_j R_j d0_j).
 *
 * Its Hessian contains mixed translation-rotation blocks from the first SO(3)
 * derivatives and local rotation-rotation blocks from the second SO(3)
 * derivatives. Only compact nodal `3 x 3` matrices are read from the active
 * workspace.
 *
 * @tparam N Number of shell nodes.
 * @param positions Current nodal midsurface positions.
 * @param rotations Compact nodal SO(3) values and derivatives.
 * @param reference_directors Nodal reference directors acted on by the rotations.
 * @param x_coefficients Interpolation coefficients of the position field.
 * @param d_coefficients Interpolation coefficients of the director field.
 * @param scale Complete stress and quadrature multiplier.
 * @param Kgeo Element geometric tangent to update.
 */
template<Index N>
void add_xd_weighted_hessian(
    const typename FRTShell<N>::MatN3& positions,
    const std::array<typename FRTShell<N>::RotationDerivatives, N>& rotations,
    const typename FRTShell<N>::MatN3& reference_directors,
    const typename FRTShell<N>::VecN&  x_coefficients,
    const typename FRTShell<N>::VecN&  d_coefficients,
    Precision                          scale,
    typename FRTShell<N>::Mat6N&       Kgeo
) {
    if (scale == Precision(0)) {
        return;
    }

    using Shell = FRTShell<N>;

    Vec3 x_value = Vec3::Zero();

    // Interpolate the current position derivative appearing in all rotational
    // blocks of this scalar product
    for (Index node = 0; node < N; ++node) {
        x_value += x_coefficients(node) * positions.row(node).transpose();
    }

    for (Index d_node = 0; d_node < N; ++d_node) {
        const Index     rot_base    = Shell::dofs_per_node * d_node + 3;
        const Precision coefficient = d_coefficients(d_node);

        if (coefficient == Precision(0)) {
            continue;
        }

        const Vec3 d0 = reference_directors.row(d_node).transpose();

        // First SO(3) derivatives generate the symmetric mixed
        // translation-rotation Hessian blocks
        for (Index a = 0; a < 3; ++a) {
            const Vec3 director_first = rotations[d_node].d1[a] * d0;

            for (Index x_node = 0; x_node < N; ++x_node) {
                const Precision x_coefficient = x_coefficients(x_node);

                if (x_coefficient == Precision(0)) {
                    continue;
                }

                const Index x_base = Shell::dofs_per_node * x_node;
                const Vec3 mixed = scale
                                 * x_coefficient
                                 * coefficient
                                 * director_first;

                Kgeo.template block<3, 1>(x_base, rot_base + a) += mixed;
                Kgeo.template block<1, 3>(rot_base + a, x_base) += mixed.transpose();
            }
        }

        // Second SO(3) derivatives contribute only to the local three-by-three
        // rotational block of the same director node
        for (Index a = 0; a < 3; ++a) {
            for (Index b = 0; b < 3; ++b) {
                const Vec3 director_second = rotations[d_node].d2[a][b] * d0;
                Kgeo(rot_base + a, rot_base + b) +=
                    scale * coefficient * x_value.dot(director_second);
            }
        }
    }
}

} // namespace

/**
 * Loads previously stored generalized shell resultants at the integration
 * points.
 *
 * The field components follow the ordering
 * `[N11,N22,N12,M11,M22,M12,Q13,Q23]`.
 *
 * @param data Active thread-local evaluation view.
 * @param ip_stress Source integration-point resultant field.
 * @param ip_start_idx First source row belonging to this element.
 */
template<Index N>
void FRTShell<N>::load_ip_resultants(
    EvaluationData& data,
    const Field&    ip_stress,
    int             ip_start_idx
) const {
    logging::error(ip_stress.components >= num_strains,
                   "FRTShell: shell resultant field requires eight components "
                   "[N11,N22,N12,M11,M22,M12,Q13,Q23]");

    // Copy the element resultants into the active reusable integration buffer
    for (Index ip = 0; ip < static_cast<Index>(data.ip_resultants.size()); ++ip) {
        const Index row = static_cast<Index>(ip_start_idx) + ip;

        for (Index component = 0; component < num_strains; ++component) {
            data.ip_resultants[static_cast<std::size_t>(ip)](component) =
                ip_stress(row, component);
        }
    }
}

/**
 * Evaluates the nonlinear shell section at every integration point.
 *
 * Generalized strains are supplied in the pointwise orthonormal reference
 * basis. The shell section receives the physical reference position and the
 * identical global basis used by the kinematic strain transformation. The
 * current consistent tangent is stored only when the caller requested B
 * matrices and therefore can consume a material tangent.
 *
 * @param data Active thread-local evaluation view.
 */
template<Index N>
void FRTShell<N>::compute_material_resultants(EvaluationData& data) const {
    logging::error(data.with_strain,
                   "FRTShell: material resultants require strain evaluation");

    ShellSection*   section = shell_section();
    const Precision scale   = topology_stiffness_scale();
    const auto&     points  = reference_data().ip_points;

    for (Index ip = 0; ip < static_cast<Index>(points.size()); ++ip) {
        const std::size_t id = static_cast<std::size_t>(ip);
        const ReferencePoint& point = points[id];
        const Vec8& strain_values = data.ip_strain[id];

        ShellGeneralizedStrain strain(strain_values);
        ShellStressResultants  resultants;
        Mat8                   tangent;

        // Construct the global pointwise material basis once for the section
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

        data.ip_resultants[id] = scale * resultants.values();

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

    const VecN dshape_r = point.dshape_rs.col(0);
    const VecN dshape_s = point.dshape_rs.col(1);
    const VecN shape    = point.shape;

    // Membrane metric Hessians contain only direct translational identities
    add_xx_weighted_hessian<N>(
        dshape_r,
        dshape_r,
        Precision(0.5) * weights(epsilon_rr),
        Kgeo
    );
    add_xx_weighted_hessian<N>(
        dshape_s,
        dshape_s,
        Precision(0.5) * weights(epsilon_ss),
        Kgeo
    );
    add_xx_weighted_hessian<N>(
        dshape_r,
        dshape_s,
        weights(gamma_rs),
        Kgeo
    );

    // Curvature Hessians contain mixed translation-rotation and local
    // rotation-rotation SO(3) blocks
    add_xd_weighted_hessian<N>(
        data.state.x,
        rotations,
        ref.d0,
        dshape_r,
        dshape_r,
        weights(kappa_rr),
        Kgeo
    );
    add_xd_weighted_hessian<N>(
        data.state.x,
        rotations,
        ref.d0,
        dshape_s,
        dshape_s,
        weights(kappa_ss),
        Kgeo
    );
    add_xd_weighted_hessian<N>(
        data.state.x,
        rotations,
        ref.d0,
        dshape_r,
        dshape_s,
        weights(kappa_rs),
        Kgeo
    );
    add_xd_weighted_hessian<N>(
        data.state.x,
        rotations,
        ref.d0,
        dshape_s,
        dshape_r,
        weights(kappa_rs),
        Kgeo
    );

    // Transverse-shear Hessians use the interpolated current director field
    add_xd_weighted_hessian<N>(
        data.state.x,
        rotations,
        ref.d0,
        dshape_r,
        shape,
        weights(gamma_r3),
        Kgeo
    );
    add_xd_weighted_hessian<N>(
        data.state.x,
        rotations,
        ref.d0,
        dshape_s,
        shape,
        weights(gamma_s3),
        Kgeo
    );
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
        const Precision t00 = points[id].invA(0, 0);
        const Precision t01 = points[id].invA(0, 1);
        const Precision t10 = points[id].invA(1, 0);
        const Precision t11 = points[id].invA(1, 1);

        StaticMatrix<3, 3> in_plane;
        in_plane << t00 * t00,               t01 * t01,               t00 * t01,
                    t10 * t10,               t11 * t11,               t10 * t11,
                    Precision(2) * t00 * t10, Precision(2) * t01 * t11,
                    t00 * t11 + t01 * t10;

        Vec8 natural_resultants = Vec8::Zero();
        natural_resultants.template segment<3>(0) =
            in_plane.transpose() * local_resultants.template segment<3>(0);
        natural_resultants.template segment<3>(3) =
            in_plane.transpose() * local_resultants.template segment<3>(3);
        natural_resultants.template segment<2>(6) =
            points[id].invA.transpose() * local_resultants.template segment<2>(6);

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
            x_a += point.dshape_a(node) * x_i;
            x_b += point.dshape_b(node) * x_i;
            a1  += point.shape(node) * rotations[node].value * point.e1;
            a2  += point.shape(node) * rotations[node].value * point.e2;
        }

        const Precision gamma_d = Precision(0.5)
                                * (a1.dot(x_b) - a2.dot(x_a));

        Vec6N B_d = Vec6N::Zero();

        // Translational derivative:
        // d(gamma_d)/du_i = 1/2 (N_i,b a1 - N_i,a a2)
        for (Index node = 0; node < num_nodes; ++node) {
            const Index base = dofs_per_node * node;
            const Vec3 derivative = Precision(0.5)
                                  * (point.dshape_b(node) * a1
                                   - point.dshape_a(node) * a2);
            B_d.template segment<3>(base) = derivative;
        }

        // Rotational derivative:
        // d(gamma_d)/dtheta_ia = N_i/2 [dR_ia e1 . x_b - dR_ia e2 . x_a]
        for (Index node = 0; node < num_nodes; ++node) {
            const Index rot_base = dofs_per_node * node + 3;
            const Precision shape = point.shape(node);

            for (Index a = 0; a < 3; ++a) {
                const Vec3 da1 = shape * rotations[node].d1[a] * point.e1;
                const Vec3 da2 = shape * rotations[node].d1[a] * point.e2;
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
        stiffness_matrix->noalias() +=
            weighted_stiffness * B_d * B_d.transpose();

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
                const Vec3 da1 = shape * rotations[rot_node].d1[a] * point.e1;
                const Vec3 da2 = shape * rotations[rot_node].d1[a] * point.e2;

                for (Index x_node = 0; x_node < num_nodes; ++x_node) {
                    const Index x_base = dofs_per_node * x_node;
                    const Vec3 mixed = Precision(0.5)
                                     * (point.dshape_b(x_node) * da1
                                      - point.dshape_a(x_node) * da2);

                    stiffness_matrix->template block<3, 1>(
                        x_base,
                        rot_base + a
                    ) += geometric_scale * mixed;
                    stiffness_matrix->template block<1, 3>(
                        rot_base + a,
                        x_base
                    ) += geometric_scale * mixed.transpose();
                }
            }

            // Pure rotation-rotation terms remain local to one nodal SO(3)
            // parameterization because the interpolated rotation field is a
            // linear sum of independent nodal rotation matrices
            for (Index a = 0; a < 3; ++a) {
                for (Index b = 0; b < 3; ++b) {
                    const Vec3 d2a1 = shape * rotations[rot_node].d2[a][b] * point.e1;
                    const Vec3 d2a2 = shape * rotations[rot_node].d2[a][b] * point.e2;
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
 * Returns the material stiffness in the active element configuration.
 *
 * The physical constitutive tangent is assembled from the current shell-section
 * response. The objective drilling tangent is added from the same kinematic
 * potential used by the nonlinear residual. Physical stress-dependent shell
 * geometric terms are intentionally omitted by this interface.
 *
 * @param buffer Caller-provided dense element matrix storage.
 * @return Mapped material and drilling tangent matrix.
 */
template<Index N>
MapMatrix FRTShell<N>::stiffness(Precision* buffer) {
    const CurrentState state = current_state();
    const EvaluationData data = init_evaluation(
        state,
        true,
        true,
        true,
        true
    );

    Mat6N Kmat;
    assemble_material_stiffness(data, Kmat);
    assemble_drill_stabilization(data, &Kmat, nullptr);

    // Remove only numerical asymmetry from the analytically symmetric tangent
    Kmat = Precision(0.5) * (Kmat + Kmat.transpose());

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = Kmat;
    return mapped;
}

/**
 * Returns the physical geometric stiffness generated by supplied shell
 * resultants.
 *
 * Drilling stabilization is not part of the physical stored resultant field
 * and is therefore intentionally excluded from this interface.
 *
 * @param buffer Caller-provided dense element matrix storage.
 * @param ip_stress Stored generalized shell resultants.
 * @param ip_start_idx First resultant row belonging to this element.
 * @return Mapped physical geometric tangent matrix.
 */
template<Index N>
MapMatrix FRTShell<N>::stiffness_geom(
    Precision*   buffer,
    const Field& ip_stress,
    int          ip_start_idx
) {
    const CurrentState state = current_state();
    const EvaluationData data = init_evaluation(
        state,
        false,
        false,
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
 * Assembles the complete consistent nonlinear tangent and matching internal
 * force.
 *
 * The supplied displacement field defines the trial configuration. The shell
 * section is updated from the current generalized strains, the physical
 * material and geometric tangents are assembled, and the objective drilling
 * force and tangent are added from the same quadratic potential. The freshly
 * evaluated generalized resultants are copied into `ip_stress_state` before the
 * element force is scattered.
 *
 * @param buffer Caller-provided dense tangent storage.
 * @param ip_stress_state Global integration-point resultant state to update.
 * @param nodal_forces Global nodal internal force field to increment.
 * @param displacement Trial nodal displacement field.
 * @return Mapped complete consistent element tangent.
 */
template<Index N>
MapMatrix FRTShell<N>::stiffness_tangent(
    Precision*   buffer,
    Field&       ip_stress_state,
    NodeData&    nodal_forces,
    const Field& displacement
) {
    logging::error(ip_stress_state.components >= num_strains,
                   "FRTShell: nonlinear stress state requires eight components");

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

    assemble_material_stiffness(data, Kmat);
    assemble_geometric_stiffness(data, Kgeo);
    assemble_internal_force(data, internal_force);

    Mat6N tangent = Kmat + Kgeo;
    assemble_drill_stabilization(data, &tangent, &internal_force);

    // Store the current generalized resultants for subsequent solver paths
    for (Index ip = 0; ip < static_cast<Index>(data.ip_resultants.size()); ++ip) {
        const Index row = this->ip_index(ip);

        for (Index component = 0; component < num_strains; ++component) {
            ip_stress_state(row, component) =
                data.ip_resultants[static_cast<std::size_t>(ip)](component);
        }
    }

    // Scatter the complete element internal force into the global nodal field
    for (Index node = 0; node < num_nodes; ++node) {
        const Index node_id = static_cast<Index>(this->node_ids[node]);

        for (Index dof = 0; dof < dofs_per_node; ++dof) {
            nodal_forces(node_id, dof) += internal_force(dofs_per_node * node + dof);
        }
    }

    tangent = Precision(0.5) * (tangent + tangent.transpose());

    MapMatrix mapped(buffer, num_dofs, num_dofs);
    mapped = tangent;
    return mapped;
}

/**
 * Recovers and scatters the nonlinear internal force from a stored generalized
 * resultant field.
 *
 * The physical shell contribution integrates `B^T n`. The objective drilling
 * force is recomputed from the current nodal configuration because it is not a
 * physical shell resultant and is intentionally not stored in the integration-
 * point stress field.
 *
 * @param node_forces Global nodal internal force field to increment.
 * @param ip_stress Stored generalized shell resultant field.
 */
template<Index N>
void FRTShell<N>::compute_internal_force_nonlinear(
    Field&       node_forces,
    const Field& ip_stress
) {
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
    assemble_drill_stabilization(data, nullptr, &internal_force);

    // Scatter the complete element force in nodal six-component ordering
    for (Index node = 0; node < num_nodes; ++node) {
        const Index node_id = static_cast<Index>(this->node_ids[node]);

        for (Index dof = 0; dof < dofs_per_node; ++dof) {
            node_forces(node_id, dof) += internal_force(dofs_per_node * node + dof);
        }
    }
}

} // namespace fem::model
