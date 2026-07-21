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

namespace {

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

template<Index N>
void add_xd_weighted_hessian(
    const std::array<typename FRTShell<N>::VectorDerivatives, N>& x_nodes,
    const std::array<typename FRTShell<N>::VectorDerivatives, N>& d_nodes,
    const typename FRTShell<N>::VecN&                            x_coefficients,
    const typename FRTShell<N>::VecN&                            d_coefficients,
    Precision                                                    scale,
    typename FRTShell<N>::Mat6N&                                 Kgeo
) {
    if (scale == Precision(0)) {
        return;
    }

    using Shell = FRTShell<N>;

    Vec3 x_value = Vec3::Zero();

    for (Index node = 0; node < N; ++node) {
        x_value += x_coefficients(node) * x_nodes[node].value;
    }

    for (Index d_node = 0; d_node < N; ++d_node) {
        const Index     rot_base      = Shell::dofs_per_node * d_node + 3;
        const Precision d_coefficient = d_coefficients(d_node);

        if (d_coefficient == Precision(0)) {
            continue;
        }

        for (Index a = 0; a < 3; ++a) {
            const Vec3 director_first = d_nodes[d_node].d1.row(rot_base + a).transpose();

            for (Index x_node = 0; x_node < N; ++x_node) {
                const Precision x_coefficient = x_coefficients(x_node);

                if (x_coefficient == Precision(0)) {
                    continue;
                }

                const Index x_base = Shell::dofs_per_node * x_node;
                const Vec3  mixed  = scale
                                  * x_coefficient
                                  * d_coefficient
                                  * director_first;

                Kgeo.template block<3, 1>(x_base, rot_base + a) += mixed;
                Kgeo.template block<1, 3>(rot_base + a, x_base) += mixed.transpose();
            }
        }

        for (Index a = 0; a < 3; ++a) {
            for (Index b = 0; b < 3; ++b) {
                Precision second = Precision(0);

                for (Index c = 0; c < 3; ++c) {
                    second += x_value(c)
                            * d_nodes[d_node].d2[c](rot_base + a, rot_base + b);
                }

                Kgeo(rot_base + a, rot_base + b) += scale * d_coefficient * second;
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
 * Pulls local generalized resultants back through the local-basis
 * transformation and the concrete MITC interpolation.
 *
 * The operation is the transpose of the local-basis transformation followed by
 * the concrete element's transposed MITC interpolation. This avoids both dense
 * Hessian fields and repeated forward probing of the MITC operator.
 */
template<Index N>
typename FRTShell<N>::Vec8 FRTShell<N>::pull_back_local_resultants(
    const ReferencePoint& point,
    const Vec8&           local_resultants
) const {
    const Precision t00 = point.invA(0, 0);
    const Precision t01 = point.invA(0, 1);
    const Precision t10 = point.invA(1, 0);
    const Precision t11 = point.invA(1, 1);

    StaticMatrix<3, 3> in_plane;
    in_plane << t00 * t00,              t01 * t01,              t00 * t01,
                t10 * t10,              t11 * t11,              t10 * t11,
                Precision(2) * t00*t10, Precision(2) * t01*t11, t00*t11 + t01*t10;

    Vec8 natural_resultants = Vec8::Zero();
    natural_resultants.template segment<3>(0) = in_plane.transpose()
        * local_resultants.template segment<3>(0);
    natural_resultants.template segment<3>(3) = in_plane.transpose()
        * local_resultants.template segment<3>(3);
    natural_resultants.template segment<2>(6) = point.invA.transpose()
        * local_resultants.template segment<2>(6);

    return natural_resultants;
}

template<Index N>
typename FRTShell<N>::GeometricStrainWeights FRTShell<N>::geometric_strain_weights(
    const ReferencePoint& point,
    const Vec8&           local_resultants
) const {
    const auto& ref = reference_data();

    GeometricStrainWeights weights;
    weights.tying.assign(ref.tying_points.size(), Vec8::Zero());

    const Vec8 natural_resultants = pull_back_local_resultants(point, local_resultants);
    pull_back_mitc_resultants(point, natural_resultants, weights.compatible, weights.tying);

    return weights;
}

/**
 * Adds a stress-weighted compatible natural Hessian directly to Kgeo.
 *
 * This contracts the eight natural strain Hessians before materializing them.
 * Only the structured translational and director-rotation blocks are touched.
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

    const VecN dshape_r = point.dshape_rs.col(0);
    const VecN dshape_s = point.dshape_rs.col(1);
    const VecN shape    = point.shape;

    add_xx_weighted_hessian<N>(dshape_r, dshape_r, Precision(0.5) * weights(epsilon_rr), Kgeo);
    add_xx_weighted_hessian<N>(dshape_s, dshape_s, Precision(0.5) * weights(epsilon_ss), Kgeo);
    add_xx_weighted_hessian<N>(dshape_r, dshape_s,                weights(gamma_rs),   Kgeo);

    add_xd_weighted_hessian<N>(data.x_nodes, data.d_nodes, dshape_r, dshape_r, weights(kappa_rr), Kgeo);
    add_xd_weighted_hessian<N>(data.x_nodes, data.d_nodes, dshape_s, dshape_s, weights(kappa_ss), Kgeo);
    add_xd_weighted_hessian<N>(data.x_nodes, data.d_nodes, dshape_r, dshape_s, weights(kappa_rs), Kgeo);
    add_xd_weighted_hessian<N>(data.x_nodes, data.d_nodes, dshape_s, dshape_r, weights(kappa_rs), Kgeo);

    add_xd_weighted_hessian<N>(data.x_nodes, data.d_nodes, dshape_r, shape, weights(gamma_r3), Kgeo);
    add_xd_weighted_hessian<N>(data.x_nodes, data.d_nodes, dshape_s, shape, weights(gamma_s3), Kgeo);
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
                   "FRTShell: geometric stiffness requires second director derivatives");
    logging::error(data.with_resultants,
                   "FRTShell: geometric stiffness requires shell resultants");

    Kgeo.setZero();

    const auto& points = reference_data().ip_points;
    const auto& tying  = reference_data().tying_points;

    for (Index ip = 0; ip < static_cast<Index>(points.size()); ++ip) {
        const std::size_t id = static_cast<std::size_t>(ip);
        const Vec8 local_resultants = data.ip_weight[id] * data.ip_resultants[id];
        const GeometricStrainWeights weights =
            geometric_strain_weights(points[id], local_resultants);

        add_weighted_natural_hessian(data, points[id], weights.compatible, Kgeo);

        for (Index tying_id = 0; tying_id < static_cast<Index>(tying.size()); ++tying_id) {
            add_weighted_natural_hessian(
                data,
                tying[static_cast<std::size_t>(tying_id)],
                weights.tying[static_cast<std::size_t>(tying_id)],
                Kgeo
            );
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
