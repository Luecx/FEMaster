/**
 * @file frt_shell_derivatives.inl
 * @brief Implements compact first and second derivatives of the fundamental
 * finite-rotation shell kinematic products.
 *
 * The finite-rotation shell strains are assembled from two recurring scalar
 * product types:
 *
 *     xx product:
 *
 *         f_xx = scale * dot(x_a, x_b),
 *
 *         x_a = sum_i a_i x_i,
 *         x_b = sum_i b_i x_i,
 *
 *     xd product:
 *
 *         f_xd = scale * dot(x_a, d_b),
 *
 *         x_a = sum_i a_i x_i,
 *         d_b = sum_i b_i d_i,
 *         d_i = R_i d0_i.
 *
 * The xx product is used for the current midsurface metric and therefore for
 * the membrane strains. Its first derivative contains only translational
 * entries, and its Hessian consists entirely of constant translational identity
 * blocks.
 *
 * The xd product is used for curvature and transverse-shear strains. Its first
 * derivative contains translational contributions from the current position
 * field and rotational contributions from the current director field. Its
 * Hessian contains symmetric translation-rotation blocks from first SO(3)
 * derivatives and local rotation-rotation blocks from second SO(3)
 * derivatives.
 *
 * All derivatives are inserted directly into the final element B matrix or
 * geometric tangent matrix. No element-sized position derivative matrices,
 * generalized strain Hessians or explicit translational identity matrices are
 * constructed or stored.
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
 * Evaluates one scaled xx product and optionally inserts its first derivative
 * directly into one row of the generalized B matrix.
 *
 * Two global vectors are formed as linear combinations of the supplied nodal
 * vector values:
 *
 *     x_a = sum_i a_i x_i,
 *     x_b = sum_i b_i x_i.
 *
 * The function returns the scalar quantity
 *
 *     f_xx = scale * dot(x_a, x_b).
 *
 * Because every nodal vector x_i depends directly and linearly on its three
 * translational degrees of freedom, the derivative with respect to x_i is
 *
 *     d f_xx / d x_i
 *         = scale * (a_i x_b + b_i x_a).
 *
 * The derivative is a three-component row vector. For all N nodes, these
 * derivatives may be represented as the N x 3 matrix
 *
 *     D_xx =
 *         scale * (a x_b^T + b x_a^T),
 *
 * whose row i contains
 *
 *     D_xx(i,:) = d f_xx / d x_i.
 *
 * When `B` is supplied, the three columns of `D_xx` are inserted into the
 * translational entries of row `B_row`. The element degrees of freedom are
 * ordered node-wise as
 *
 *     q_i =
 *     [u_x,i, u_y,i, u_z,i, theta_x,i, theta_y,i, theta_z,i].
 *
 * Consequently, the translational components of all nodes occupy the strided
 * column sequences
 *
 *     u_x: 0, 6, 12, ...
 *     u_y: 1, 7, 13, ...
 *     u_z: 2, 8, 14, ...
 *
 * For example, for a four-node shell, the selected B row receives
 *
 *     [df/du_x,0  df/du_y,0  df/du_z,0  0  0  0
 *      df/du_x,1  df/du_y,1  df/du_z,1  0  0  0
 *      df/du_x,2  df/du_y,2  df/du_z,2  0  0  0
 *      df/du_x,3  df/du_y,3  df/du_z,3  0  0  0].
 *
 * The rotational entries are not modified because the xx product depends only
 * on the supplied nodal vector values.
 *
 * @tparam N Number of shell nodes.
 * @param nodal_values Nodal global vectors x_i.
 * @param a_coefficients Coefficients a_i forming x_a.
 * @param b_coefficients Coefficients b_i forming x_b.
 * @param scale Scalar multiplier applied to the value and its derivative.
 * @param B_row Generalized strain row receiving the first derivative.
 * @param B Optional generalized B matrix to update.
 * @return Scaled scalar product `scale * dot(x_a, x_b)`.
 */
template<Index N>
Precision evaluate_xx_product(
    const typename FRTShell<N>::MatN3& nodal_values,
    const typename FRTShell<N>::VecN&  a_coefficients,
    const typename FRTShell<N>::VecN&  b_coefficients,
    Precision                          scale,
    Index                              B_row,
    typename FRTShell<N>::Mat8x6N*     B
) {
    using Shell = FRTShell<N>;

    // Form the two interpolated global vectors:
    //
    //     x_a = sum_i a_i x_i,
    //     x_b = sum_i b_i x_i.
    const Vec3 x_a = nodal_values.transpose() * a_coefficients;
    const Vec3 x_b = nodal_values.transpose() * b_coefficients;

    if (B) {
        // Assemble the derivatives with respect to all nodal vectors at once:
        //
        //     D_xx(i,:) = scale * (a_i x_b + b_i x_a)^T.
        //
        // The two outer products produce the complete N x 3 derivative matrix
        // directly, without a loop over the shell nodes.
        const typename Shell::MatN3 nodal_derivative =
            scale * (a_coefficients * x_b.transpose() + b_coefficients * x_a.transpose());

        // Scatter the first derivative component of every node into the global
        // x-translation columns 0, 6, 12, ... of the selected B row.
        B->row(B_row)(Eigen::seqN(0, N, Shell::dofs_per_node)) += nodal_derivative.col(0).transpose();

        // Scatter the second derivative component of every node into the global
        // y-translation columns 1, 7, 13, ... of the selected B row.
        B->row(B_row)(Eigen::seqN(1, N, Shell::dofs_per_node)) += nodal_derivative.col(1).transpose();

        // Scatter the third derivative component of every node into the global
        // z-translation columns 2, 8, 14, ... of the selected B row.
        B->row(B_row)(Eigen::seqN(2, N, Shell::dofs_per_node)) += nodal_derivative.col(2).transpose();
    }

    return scale * x_a.dot(x_b);
}

/**
 * Evaluates one scaled xd product and optionally inserts its first derivative
 * directly into one row of the generalized B matrix.
 *
 * The current position vector is formed from the nodal midsurface positions:
 *
 *     x_a = sum_i a_i x_i.
 *
 * The current director vector is formed from the rotated nodal reference
 * directors:
 *
 *     d_b = sum_i b_i d_i,
 *
 * with
 *
 *     d_i = R_i d0_i.
 *
 * The function returns
 *
 *     f_xd = scale * dot(x_a, d_b).
 *
 * The first derivative consists of two separate parts.
 *
 * The derivative with respect to the translational vector x_i follows directly
 * from the linear interpolation of x_a:
 *
 *     d f_xd / d x_i
 *         = scale * a_i d_b.
 *
 * For all nodes simultaneously, the N x 3 translational derivative matrix is
 *
 *     D_x =
 *         scale * a d_b^T.
 *
 * Its three columns are inserted into the strided translational entries of the
 * selected B row.
 *
 * The director at node i depends on the three local axis-angle coordinates
 * theta_i,c through
 *
 *     d_i = R_i(theta_i) d0_i.
 *
 * Its first rotational derivative is
 *
 *     d d_i / d theta_i,c
 *         = (dR_i / d theta_i,c) d0_i.
 *
 * Therefore,
 *
 *     d f_xd / d theta_i,c
 *         = scale
 *         * b_i
 *         * dot(x_a, (dR_i / d theta_i,c) d0_i).
 *
 * These values are inserted into the three rotational entries belonging to the
 * same node:
 *
 *     6i + 3,
 *     6i + 4,
 *     6i + 5.
 *
 * Unlike the translational derivative, the rotational derivative cannot be
 * assembled through one common outer product. Every node and every rotational
 * component uses a different SO(3) derivative matrix.
 *
 * @tparam N Number of shell nodes.
 * @param data Active evaluation data containing current positions, directors
 * and optional first SO(3) derivatives.
 * @param reference_directors Nodal reference directors d0_i.
 * @param x_coefficients Coefficients a_i forming x_a.
 * @param d_coefficients Coefficients b_i forming d_b.
 * @param scale Scalar multiplier applied to the value and its derivative.
 * @param B_row Generalized strain row receiving the first derivative.
 * @param B Optional generalized B matrix to update.
 * @return Scaled scalar product `scale * dot(x_a, d_b)`.
 */
template<Index N>
Precision evaluate_xd_product(
    const typename FRTShell<N>::EvaluationData& data,
    const typename FRTShell<N>::MatN3&          reference_directors,
    const typename FRTShell<N>::VecN&           x_coefficients,
    const typename FRTShell<N>::VecN&           d_coefficients,
    Precision                                   scale,
    Index                                       B_row,
    typename FRTShell<N>::Mat8x6N*              B
) {
    using Shell = FRTShell<N>;

    // Form the current position and director vectors:
    //
    //     x_a = sum_i a_i x_i,
    //     d_b = sum_i b_i d_i.
    const Vec3 x_a = data.state.x.transpose() * x_coefficients;
    const Vec3 d_b = data.state.d.transpose() * d_coefficients;

    if (B) {
        logging::error(data.rotations != nullptr,
                       "FRTShell: B evaluation requires nodal rotation derivatives");

        // Assemble all translational derivatives simultaneously:
        //
        //     D_x(i,:) = scale * a_i d_b^T.
        //
        // The outer product contains one derivative row for every nodal
        // translational vector.
        const typename Shell::MatN3 translational_derivative =
            scale * x_coefficients * d_b.transpose();

        // Scatter the first derivative component of every node into the global
        // x-translation columns 0, 6, 12, ... of the selected B row.
        B->row(B_row)(Eigen::seqN(0, N, Shell::dofs_per_node)) += translational_derivative.col(0).transpose();

        // Scatter the second derivative component of every node into the global
        // y-translation columns 1, 7, 13, ... of the selected B row.
        B->row(B_row)(Eigen::seqN(1, N, Shell::dofs_per_node)) += translational_derivative.col(1).transpose();

        // Scatter the third derivative component of every node into the global
        // z-translation columns 2, 8, 14, ... of the selected B row.
        B->row(B_row)(Eigen::seqN(2, N, Shell::dofs_per_node)) += translational_derivative.col(2).transpose();

        const auto& rotations = *data.rotations;

        // Evaluate the three rotational derivatives independently for every
        // nodal director. Each director depends only on the rotational
        // coordinates of its own node, so no off-node rotational entries occur
        // in the first derivative.
        for (Index node = 0; node < N; ++node) {
            const Precision coefficient = scale * d_coefficients(node);

            if (coefficient == Precision(0)) {
                continue;
            }

            const Index rot_base = Shell::dofs_per_node * node + 3;
            const Vec3  d0       = reference_directors.row(node).transpose();

            for (Index component = 0; component < 3; ++component) {
                // Differentiate the current nodal director:
                //
                //     d d_i / d theta_i,c
                //         = (dR_i / d theta_i,c) d0_i.
                const Vec3 director_derivative = rotations[node].d1[component] * d0;

                // Insert
                //
                //     scale * b_i * dot(x_a, d d_i / d theta_i,c)
                //
                // into the corresponding rotational entry of the selected
                // generalized B row.
                (*B)(B_row, rot_base + component) += coefficient * x_a.dot(director_derivative);
            }
        }
    }

    return scale * x_a.dot(d_b);
}

/**
 * Adds the scaled Hessian of one xx product directly to the element tangent
 * matrix.
 *
 * The underlying scalar quantity is
 *
 *     f_xx = scale * dot(x_a, x_b),
 *
 * with
 *
 *     x_a = sum_i a_i x_i,
 *     x_b = sum_i b_i x_i.
 *
 * Its first derivative with respect to nodal vector x_i is
 *
 *     d f_xx / d x_i
 *         = scale * (a_i x_b + b_i x_a).
 *
 * Differentiating again with respect to x_j gives
 *
 *     d² f_xx / (d x_i d x_j)
 *         = scale * (a_i b_j + b_i a_j) I.
 *
 * Every nodal pair therefore contributes one 3 x 3 translational identity
 * block. The block is multiplied by the scalar
 *
 *     block_scale
 *         = scale * (a_i b_j + b_i a_j).
 *
 * The element degrees of freedom are ordered as
 *
 *     [u_x,i, u_y,i, u_z,i, theta_x,i, theta_y,i, theta_z,i]
 *
 * for every node i. Consequently, the Hessian modifies only the upper-left
 * translational 3 x 3 part of each nodal pair block:
 *
 *            node j
 *          translation   rotation
 *
 * node i   [ c I           0    ]
 *          [  0            0    ],
 *
 * where
 *
 *     c = scale * (a_i b_j + b_i a_j).
 *
 * No explicit identity matrix and no complete intermediate Hessian are
 * constructed. Adding the scalar to the diagonal of the destination 3 x 3
 * block produces the required identity contribution directly.
 *
 * The double loop over nodal pairs is intrinsic to the destination layout:
 * every pair `(i,j)` occupies a different strided 3 x 3 block in the complete
 * 6N x 6N element matrix.
 *
 * @tparam N Number of shell nodes.
 * @param a_coefficients Coefficients a_i forming x_a.
 * @param b_coefficients Coefficients b_i forming x_b.
 * @param scale Complete scalar multiplier, typically containing a generalized
 * resultant weight and the reference-area quadrature factor.
 * @param Kgeo Element tangent matrix receiving the scaled Hessian.
 */
template<Index N>
void add_xx_hessian(
    const typename FRTShell<N>::VecN& a_coefficients,
    const typename FRTShell<N>::VecN& b_coefficients,
    Precision                         scale,
    typename FRTShell<N>::Mat6N&      Kgeo
) {
    using Shell = FRTShell<N>;

    if (scale == Precision(0)) {
        return;
    }

    // Insert the constant translational identity block associated with every
    // ordered pair of nodal vectors x_i and x_j.
    for (Index i = 0; i < N; ++i) {
        const Index row = Shell::dofs_per_node * i;

        for (Index j = 0; j < N; ++j) {
            const Index col = Shell::dofs_per_node * j;

            // Scalar multiplying the 3 x 3 identity block:
            //
            //     scale * (a_i b_j + b_i a_j).
            const Precision block_scale =
                scale * (a_coefficients(i) * b_coefficients(j) + b_coefficients(i) * a_coefficients(j));

            if (block_scale == Precision(0)) {
                continue;
            }

            // Add block_scale * I directly to the translational block. Only
            // its three diagonal entries are modified; the rotational rows and
            // columns remain untouched.
            Kgeo.template block<3, 3>(row, col).diagonal().array() += block_scale;
        }
    }
}

/**
 * Adds the scaled Hessian of one xd product directly to the element tangent
 * matrix.
 *
 * The underlying scalar quantity is
 *
 *     f_xd = scale * dot(x_a, d_b),
 *
 * where
 *
 *     x_a = sum_i a_i x_i,
 *
 *     d_b = sum_j b_j d_j,
 *
 *     d_j = R_j d0_j.
 *
 * Since x_a depends linearly on nodal translations and d_b depends nonlinearly
 * on the nodal rotations, the Hessian contains two types of nonzero blocks:
 *
 *  1. translation-rotation blocks generated by first SO(3) derivatives,
 *  2. rotation-rotation blocks generated by second SO(3) derivatives.
 *
 * There are no translation-translation blocks because
 *
 *     d² x_a / (d x_i d x_j) = 0,
 *
 * and there are no rotation-rotation blocks coupling two different director
 * nodes because each nodal director d_j depends only on its own nodal rotation
 * vector theta_j.
 *
 * For translation node i and director node j,
 *
 *     d² f_xd / (d x_i d theta_j,c)
 *         = scale
 *         * a_i
 *         * b_j
 *         * (dR_j / d theta_j,c) d0_j.
 *
 * This is a three-component column and is added to the block
 *
 *     Kgeo(
 *         translational rows of node i,
 *         rotational component c of node j
 *     ).
 *
 * The transpose is added to the corresponding rotation-translation block so
 * the analytically symmetric scalar Hessian is assembled directly:
 *
 *     K_xθ += mixed,
 *     K_θx += mixed^T.
 *
 * The pure rotational block of director node j follows from
 *
 *     d² d_j / (d theta_j,a d theta_j,b)
 *         = (d²R_j / d theta_j,a d theta_j,b) d0_j.
 *
 * Therefore,
 *
 *     d² f_xd / (d theta_j,a d theta_j,b)
 *         = scale
 *         * b_j
 *         * dot(
 *               x_a,
 *               (d²R_j / d theta_j,a d theta_j,b) d0_j
 *           ).
 *
 * These nine scalar entries are added only to the local 3 x 3 rotational block
 * of node j.
 *
 * The current interpolated position vector x_a is evaluated once before the
 * nodal loops because it is reused by every rotation-rotation contribution.
 * Only the compact nodal 3 x 3 SO(3) derivative matrices are accessed. No dense
 * generalized strain Hessian is materialized.
 *
 * @tparam N Number of shell nodes.
 * @param positions Current nodal midsurface positions x_i.
 * @param rotations Nodal SO(3) rotations and their first and second derivatives.
 * @param reference_directors Nodal reference directors d0_i.
 * @param x_coefficients Coefficients a_i forming x_a.
 * @param d_coefficients Coefficients b_i forming d_b.
 * @param scale Complete scalar multiplier, typically containing a generalized
 * resultant weight and the reference-area quadrature factor.
 * @param Kgeo Element tangent matrix receiving the scaled Hessian.
 */
template<Index N>
void add_xd_hessian(
    const typename FRTShell<N>::MatN3& positions,
    const std::array<typename FRTShell<N>::RotationDerivatives, N>& rotations,
    const typename FRTShell<N>::MatN3& reference_directors,
    const typename FRTShell<N>::VecN&  x_coefficients,
    const typename FRTShell<N>::VecN&  d_coefficients,
    Precision                          scale,
    typename FRTShell<N>::Mat6N&       Kgeo
) {
    using Shell = FRTShell<N>;

    if (scale == Precision(0)) {
        return;
    }

    // Interpolate the current position vector once:
    //
    //     x_a = sum_i a_i x_i.
    //
    // The result is reused by all rotation-rotation Hessian entries.
    const Vec3 x_a = positions.transpose() * x_coefficients;

    // Process every nodal director contributing to d_b. Different director
    // nodes are rotationally independent, so each node contributes one set of
    // mixed blocks and one local rotational block.
    for (Index d_node = 0; d_node < N; ++d_node) {
        const Precision d_coefficient = d_coefficients(d_node);

        if (d_coefficient == Precision(0)) {
            continue;
        }

        const Index rot_base = Shell::dofs_per_node * d_node + 3;
        const Vec3  d0       = reference_directors.row(d_node).transpose();

        // ---------------------------------------------------------------------
        // Mixed translation-rotation Hessian blocks
        // ---------------------------------------------------------------------
        //
        // The first derivative of the current director with respect to one
        // nodal rotation component is
        //
        //     d d_j / d theta_j,a
        //         = (dR_j / d theta_j,a) d0_j.
        for (Index component = 0; component < 3; ++component) {
            const Vec3 director_first = rotations[d_node].d1[component] * d0;

            // Every position node contributing to x_a couples to the current
            // director node through one symmetric translation-rotation block.
            for (Index x_node = 0; x_node < N; ++x_node) {
                const Precision x_coefficient = x_coefficients(x_node);

                if (x_coefficient == Precision(0)) {
                    continue;
                }

                const Index x_base = Shell::dofs_per_node * x_node;

                // Three-component mixed derivative:
                //
                //     scale * a_i * b_j
                //           * (dR_j / d theta_j,c) d0_j.
                const Vec3 mixed = scale * x_coefficient * d_coefficient * director_first;

                // Add the translation-rotation block and its symmetric
                // rotation-translation counterpart.
                Kgeo.template block<3, 1>(x_base              , rot_base + component) += mixed;
                Kgeo.template block<1, 3>(rot_base + component, x_base              ) += mixed.transpose();
            }
        }

        // ---------------------------------------------------------------------
        // Local rotation-rotation Hessian block
        // ---------------------------------------------------------------------
        //
        // Second SO(3) derivatives remain local to the current director node.
        // The nine combinations of nodal rotation components form its complete
        // 3 x 3 rotational Hessian block.
        for (Index a = 0; a < 3; ++a) {
            for (Index b = 0; b < 3; ++b) {
                // Second derivative of the current nodal director:
                //
                //     d² d_j / (d theta_j,a d theta_j,b)
                //         = (d²R_j / d theta_j,a d theta_j,b) d0_j.
                const Vec3 director_second = rotations[d_node].d2[a][b] * d0;

                // Add the contracted second director derivative:
                //
                //     scale * b_j * dot(x_a, director_second).
                Kgeo(rot_base + a, rot_base + b) += scale * d_coefficient * x_a.dot(director_second);
            }
        }
    }
}

} // namespace

} // namespace fem::model