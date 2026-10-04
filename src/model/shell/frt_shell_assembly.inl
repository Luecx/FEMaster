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
 * material history; base-state stiffness and perturbation geometric stiffness are state-neutral.
 *
 * @see FRTShell
 *
 * @author Finn Eggers
 * @date 21.07.2026
 */

#include "frt_shell.h"

#include "../../core/logging.h"

#include <cmath>
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

        // Apply the transpose of the pointwise natural-to-local strain map.
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

        // The forward path deskews compatible strains independently at every
        // sampling point before MITC interpolation. Apply the exact transposed
        // pointwise maps here after the MITC pull-back so the raw compatible
        // strain Hessians receive work-conjugate weights.
        compatible_weights =
            director_deskew_natural_transform(points[id]).transpose()
            * compatible_weights;

        // Add the compatible integration-point Hessian contribution
        add_weighted_natural_hessian(data, points[id], compatible_weights, Kgeo);

        // Add all compatible tying-point Hessian contributions
        for (Index tying_id = 0; tying_id < static_cast<Index>(tying.size()); ++tying_id) {
            const std::size_t tying_index = static_cast<std::size_t>(tying_id);
            Vec8& tying_weight = data.geometric_tying_weights[tying_index];
            tying_weight =
                director_deskew_natural_transform(tying[tying_index]).transpose()
                * tying_weight;

            add_weighted_natural_hessian(
                data,
                tying[tying_index],
                tying_weight,
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
 * Maps isotropic thermal free expansion into the assumed shell strain space.
 *
 * Unit midsurface dilation u_i = X_i with unchanged nodal directors gives the
 * compatible reference metric and curvature increments. Deskewing and applying
 * the topology-specific MITC operator must match the mechanical strain path:
 * on curved MITC8 elements its assumed field differs from the pointwise field.
 * The same initial strain is used for loads, stress recovery and thermal stress increments.
 *
 * The supplied scalar retains the existing pointwise temperature interpolation.
 * This models constant-through-thickness expansion and introduces no thermal
 * director rotations or material-state updates.
 *
 * @param point Reference geometry at the target integration or recovery point.
 * @param free_strain Interpolated isotropic expansion alpha times DeltaT.
 * @return Thermal generalized strain in the local orthonormal reference basis.
 */
template<Index N>
typename FRTShell<N>::Vec8 FRTShell<N>::thermal_generalized_strain(
    const ReferencePoint& point,
    Precision             free_strain
) const {
    // Reference-temperature points have no thermal initial strain
    if (free_strain == Precision(0)) {
        return Vec8::Zero();
    }

    // Uniform reference dilation changes the metric and curvature, while the
    // nodal directors retain their orientation. The raw covariant shear change
    // is removed by the same reference-director deskew used mechanically.
    const auto compatible_dilation = [&](const ReferencePoint& sample) {
        Vec8 strain = Vec8::Zero();
        strain(0) = sample.X_rs.col(0).squaredNorm();
        strain(1) = sample.X_rs.col(1).squaredNorm();
        strain(2) = Precision(2) * sample.X_rs.col(0).dot(sample.X_rs.col(1));
        strain(3) = sample.X_rs.col(0).dot(sample.D_rs.col(0));
        strain(4) = sample.X_rs.col(1).dot(sample.D_rs.col(1));
        strain(5) = sample.X_rs.col(0).dot(sample.D_rs.col(1))
                  + sample.X_rs.col(1).dot(sample.D_rs.col(0));
        strain(6) = sample.X_rs.col(0).dot(sample.D);
        strain(7) = sample.X_rs.col(1).dot(sample.D);
        return (director_deskew_natural_transform(sample) * strain).eval();
    };

    // Keep tying values in local storage so result recovery cannot invalidate
    // the active mechanical evaluation's thread-local workspace.
    const auto& tying_points = reference_data().tying_points;
    std::vector<Vec8> tying_strain(tying_points.size());
    for (std::size_t tying = 0; tying < tying_points.size(); ++tying) {
        tying_strain[tying] = compatible_dilation(tying_points[tying]);
    }

    EvaluationData thermal_data;
    thermal_data.with_strain      = true;
    thermal_data.tying_strain_nat = Span<Vec8>(tying_strain);

    // Apply the identical assumed-strain interpolation and local basis mapping
    // before scaling by the temperature-induced expansion.
    Vec8 thermal_strain = compatible_dilation(point);
    apply_mitc_natural(thermal_data, point, thermal_strain, nullptr);
    transform_strain_to_local(point, thermal_strain, nullptr);
    return free_strain * thermal_strain;
}

/**
 * Integrates equivalent nodal forces from a scalar midsurface temperature field.
 *
 * The prescribed temperature is constant through the section thickness. For
 * isotropic thermal expansion the free generalized membrane strain is equal in
 * both tangent directions. On a curved reference surface it also contains the
 * change of the generalized curvature produced by uniform midsurface scaling:
 * epsilon_th times X_,a dot D_,b. Without this term a uniformly heated cylinder
 * can avoid artificial bending energy by opening its seam instead of expanding
 * radially. The free transverse-shear strain remains zero. Multiplying by the
 * complete section tangent (including ABD coupling) yields generalized thermal
 * membrane forces and bending moments. The reference MITC B matrix and the
 * linear-stiffness integration weights map these resultants into consistent
 * forces and moments at all six nodal DOFs. Thermal metric and curvature
 * increments pass through the same MITC operator as mechanical strains.
 *
 * This is the linear/reference thermal RHS, not a constitutive update. It does
 * not modify material state; finite-state thermal constitutive response still
 * requires separate thermal state handling.
 */
template<Index N>
void FRTShell<N>::apply_tload(Field& node_loads, const Field& node_temp, Precision ref_temp) {
    logging::error(node_temp.domain == FieldDomain::NODE && node_temp.components == 1,
                   "FRTShell: thermal loading requires a scalar nodal temperature field");
    logging::error(node_loads.domain == FieldDomain::NODE
                   && node_loads.components >= dofs_per_node,
                   "FRTShell: thermal loading requires six nodal load components");
    logging::error(std::isfinite(ref_temp),
                   "FRTShell: thermal reference temperature must be finite");

    const auto material = this->get_material();
    logging::error(material->has_thermal_expansion(),
                   "FRTShell: material has no thermal expansion for element ", this->elem_id);

    VecN nodal_temperatures;
    for (Index node = 0; node < num_nodes; ++node) {
        const Index node_id = static_cast<Index>(this->node_ids[node]);
        const Precision temperature = node_temp(node_id, 0);
        // Match the solid TLOAD convention for undefined nodal temperatures.
        nodal_temperatures(node) = std::isfinite(temperature) ? temperature : ref_temp;
    }

    const Precision alpha = material->get_thermal_expansion();
    const EvaluationData data = init_evaluation(
        reference_state(), true, true, false, false, false
    );
    const auto& points = reference_data().ip_points;
    Vec6N thermal_force = Vec6N::Zero();

    for (Index ip = 0; ip < static_cast<Index>(points.size()); ++ip) {
        const std::size_t id = static_cast<std::size_t>(ip);
        const ReferencePoint& point = points[id];
        const Precision temperature =
            shape_function(point.r, point.s).dot(nodal_temperatures);
        const Precision free_strain = alpha * (temperature - ref_temp);

        const Vec8 thermal_strain = thermal_generalized_strain(point, free_strain);

        thermal_force.noalias() += (point.w * point.detJ)
            * data.ip_B[id].transpose()
            * (data.ip_tangent[id] * thermal_strain);
    }

    for (Index node = 0; node < num_nodes; ++node) {
        const Index node_id = static_cast<Index>(this->node_ids[node]);
        for (Index dof = 0; dof < dofs_per_node; ++dof) {
            node_loads(node_id, dof) += thermal_force(dofs_per_node * node + dof);
        }
    }
}

template<Index N>
void FRTShell<N>::apply_thermal_free_strain(Field& thermal_free_strain,
                                            const Field& node_temp,
                                            Precision ref_temp) {
    logging::error(thermal_free_strain.domain == FieldDomain::ELEMENT_NODAL
                   && thermal_free_strain.components == 1,
                   "FRTShell: thermal free strain requires scalar ELEMENT_NODAL storage");
    logging::error(node_temp.domain == FieldDomain::NODE && node_temp.components == 1,
                   "FRTShell: thermal free strain requires a scalar nodal temperature field");
    logging::error(std::isfinite(ref_temp),
                   "FRTShell: thermal reference temperature must be finite");

    const auto material = this->get_material();
    logging::error(material->has_thermal_expansion(),
                   "FRTShell: material has no thermal expansion for element ", this->elem_id);

    const Precision alpha = material->get_thermal_expansion();
    for (Index node = 0; node < num_nodes; ++node) {
        const Precision value =
            node_temp(static_cast<Index>(this->node_ids[node]), 0);
        const Precision temperature = std::isfinite(value) ? value : ref_temp;
        thermal_free_strain(static_cast<Index>(this->elem_nodal_offset) + node, 0) +=
            alpha * (temperature - ref_temp);
    }
}

/**
 * Evaluates shell force and tangent operators about a selected base state u0.
 *
 * A null linearization denotes the reference state u0 = 0. A non-null
 * linearization activates the finite-rotation shell kinematics at the supplied
 * base state. The complete tangent is always evaluated at u0 and contains the
 * material and geometric contributions generated by the resultants already
 * present at that state.
 *
 * Internal force at a different requested state u is returned through the
 * first-order approximation
 *
 *     f(u) = f(u0) + K_T(u0) (u - u0).
 *
 * The separately requested geometric stiffness is generated only by the
 * linearized resultant increment from u0 to u. Existing resultants at u0 remain
 * part of K_T(u0) and are not repeated in the separate geometric operator.
 *
 * Constitutive trial state can only be written for an exact physical evaluation
 * at u0.
 */
template<Index N>
MapMatrix FRTShell<N>::evaluate(
    Precision*   tangent_buffer,
    Precision*   geometric_tangent_buffer,
    NodeData*    internal_force_output,
    const Field* displacement,
    const Field* linearization,
    const Field* thermal_free_strain,
    bool         update_state
) {
    const bool with_tangent   = tangent_buffer           != nullptr;
    const bool with_geometric = geometric_tangent_buffer != nullptr;
    const bool with_force     = internal_force_output    != nullptr;

    if (!with_tangent && !with_geometric && !with_force) {
        return MapMatrix(nullptr, 0, 0);
    }

    logging::error(!with_force || displacement != nullptr,
        "FRTShell: internal force evaluation requires displacement");
    logging::error(!with_force || internal_force_output->components >= dofs_per_node,
        "FRTShell: internal force requires six nodal components");
    logging::error(!update_state || (linearization != nullptr && displacement == linearization),
        "FRTShell: material state requires an exact evaluation at the linearization state");
    logging::error(!thermal_free_strain || linearization == nullptr,
        "FRTShell: finite-rotation thermal free strain is not implemented");
    logging::error(!thermal_free_strain || (thermal_free_strain->domain == FieldDomain::ELEMENT_NODAL && thermal_free_strain->components == 1),
        "FRTShell: thermal free strain must be scalar ELEMENT_NODAL data");

    // -------------------------------------------------------------------------
    // exact state at the linearization point u0
    // -------------------------------------------------------------------------

    // A null linearization denotes q0 = 0. Reference and finite base states use
    // the same nonlinear shell kinematics and the same Green-Lagrange section
    // response; only the supplied nodal state differs.
    const Vec6N q0 = linearization ? element_displacement_vector(*linearization) : Vec6N::Zero();
    const Vec6N q  = displacement ? element_displacement_vector(*displacement) : q0;
    const Vec6N delta = q - q0;

    const CurrentState state = linearization ? current_state_from_displacement(*linearization) : reference_state();

    // Internal force away from q0 is continued affinely:
    //
    //     f(q) ~= f(q0) + K_T(q0) (q - q0).
    const bool affine_force          = with_force && displacement != linearization;
    const bool need_complete_tangent = with_tangent || affine_force;
    const bool need_B                = with_force || need_complete_tangent || with_geometric;
    const bool need_G                = need_complete_tangent || with_geometric;
    const bool need_resultants       = with_force || need_complete_tangent || with_geometric;

    EvaluationData data = init_evaluation(state, true, need_B, need_G, need_resultants, update_state);

    Vec6N force = Vec6N::Zero();
    if (with_force) {
        assemble_internal_force(data, force);
    }

    // The complete tangent at q0 contains
    //
    //     K_T(q0) = K_M(q0) + K_G(n0),
    //
    // where n0 are the exact generalized resultants at the base state.
    Mat6N complete  = Mat6N::Zero();
    Mat6N geometric = Mat6N::Zero();

    if (need_complete_tangent) {
        Mat6N material;
        Mat6N geometric_base;

        assemble_material_stiffness(data, material);
        assemble_geometric_stiffness(data, geometric_base);
        complete = material + geometric_base;

        assemble_drill_stabilization(data, &complete, with_force ? &force : nullptr);
    } else if (with_force) {
        assemble_drill_stabilization(data, nullptr, &force);
    }

    // -------------------------------------------------------------------------
    // linear continuation from u0 to u
    // -------------------------------------------------------------------------

    // Linearization about q0 gives
    //
    //     Delta epsilon = B0 Delta q,
    //     Delta n       = H0 Delta epsilon.
    //
    // Reference thermal free strain is treated as an additional negative
    // generalized strain increment.
    Vec6N thermal_force = Vec6N::Zero();
    thread_local std::vector<Vec8> resultant_increment;

    if (with_geometric) {
        resultant_increment.resize(reference_data().ip_points.size());
    }

    if (with_geometric || (with_force && thermal_free_strain)) {
        const auto& points = reference_data().ip_points;

        for (Index ip = 0; ip < static_cast<Index>(points.size()); ++ip) {
            const std::size_t id = static_cast<std::size_t>(ip);
            const ReferencePoint& point = points[id];

            Vec8 increment = data.ip_tangent[id] * (data.ip_B[id] * delta);

            if (thermal_free_strain) {
                const Precision free_strain      = thermal_free_strain_at(thermal_free_strain, point.r, point.s);
                const Vec8 thermal_strain        = thermal_generalized_strain(point, free_strain);
                const Vec8 thermal_resultant     = data.ip_tangent[id] * thermal_strain;

                increment -= thermal_resultant;

                if (with_force) {
                    thermal_force.noalias() += (point.w * point.detJ) * data.ip_B[id].transpose() * thermal_resultant;
                }
            }

            if (with_geometric) {
                resultant_increment[id] = increment;
            }
        }
    }

    // The separate geometric operator contains only the perturbation resultants
    // Delta n. Base resultants n0 are already contained in K_T(q0).
    if (with_geometric) {
        EvaluationData geometric_data = data;
        geometric_data.ip_resultants   = Span<Vec8>(resultant_increment);
        geometric_data.with_resultants = true;

        assemble_geometric_stiffness(geometric_data, geometric);

        MapMatrix mapped(geometric_tangent_buffer, num_dofs, num_dofs);
        mapped = geometric;
    }

    if (affine_force) {
        force.noalias() += complete * delta;
    }
    if (with_force && thermal_free_strain) {
        force.noalias() -= thermal_force;
    }

    if (with_force) {
        for (Index node = 0; node < num_nodes; ++node) {
            const Index node_id = static_cast<Index>(this->node_ids[node]);
            for (Index dof = 0; dof < dofs_per_node; ++dof) {
                (*internal_force_output)(node_id, dof) += force(dofs_per_node * node + dof);
            }
        }
    }

    if (with_tangent) {
        MapMatrix mapped(tangent_buffer, num_dofs, num_dofs);
        mapped = complete;
    }

    if (with_tangent) return MapMatrix(tangent_buffer, num_dofs, num_dofs);
    if (with_geometric) return MapMatrix(geometric_tangent_buffer, num_dofs, num_dofs);
    return MapMatrix(nullptr, 0, 0);
}

} // namespace fem::model
