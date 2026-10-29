/**
 * @file eigenvalue_estimates.cpp
 * @brief Implements inexpensive eigenvalue estimates in the matrix-tools layer.
 *
 * Structural load cases supply reduced symmetric operators for
 *
 *     A q = lambda B q,
 *
 * where A = T^T K T is elastic stiffness. For modal analysis B = T^T M T and
 * lambda = omega^2; for linear buckling B = -T^T Kg T and lambda is the load
 * multiplier. The geometric-stiffness sign is supplied by the caller.
 *
 * The Ritz estimator samples quadratic global displacement fields at independent
 * master DOFs. ConstraintMap supplies their coordinate interpretation through
 * u = T q. Six global rigid-motion candidates are checked against A and used
 * to remove rigid content. Elastic-energy orthonormalization then produces a
 * small trial space for a reciprocal generalized eigenproblem. Using A as the
 * positive metric permits B to be indefinite in buckling problems.
 *
 * The separate diagonal-sum estimator supplies a geometry-independent modal
 * fallback when the quadratic space is empty or otherwise unusable. Both
 * estimates leave input matrices, constraints and model state unchanged.
 * Constraint construction, assembly, fallback selection and spectral shifts
 * remain responsibilities of the calling load case; this file performs no
 * factorization of the full reduced system.
 *
 * @see eigenvalue_estimates.h
 * @see constraint::ConstraintMap
 *
 * @author Finn Eggers
 * @date 29.10.2026
 */

#include "eigenvalue_estimates.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace fem {
namespace mattools {

/**
 * @brief Computes a modal eigenvalue scale without geometric trial fields.
 *
 * For positive semidefinite stiffness A and mass B, use
 *
 *     lambda_scale = trace(A) / trace(B) = sum_i A_ii / sum_i B_ii.
 *
 * This has the units of omega^2 and is positive whenever the operators contain
 * elastic stiffness and mass. The sums are formed before division, so a DOF
 * with B_ii = 0 does not cause an individual division by zero. Coordinates,
 * polynomial samples and the stiffness nullspace are not needed.
 *
 * This ratio measures a characteristic scale in the supplied coordinates. It
 * neither bounds nor specifically estimates the smallest flexible eigenvalue,
 * and it is not invariant under arbitrary rescaling of individual coordinates.
 * The caller decides whether it is appropriate as a fallback and chooses its
 * shift. An indefinite buckling B is outside this estimator's intended use;
 * positive diagonal sums alone would not establish positive semidefiniteness.
 *
 * @param A Reduced symmetric positive semidefinite elastic stiffness.
 * @param B Reduced symmetric positive semidefinite mass in the same coordinates.
 * @return Finite positive diagonal-sum ratio, or NaN for incompatible dimensions,
 *         nonfinite arithmetic or nonpositive stiffness/mass sums. No input is
 *         modified and no factorization is performed.
 */
Precision estimate_spectral_scale(const SparseMatrix& A, const SparseMatrix& B) {
    // ---------------------------------------------------------------------
    // (1) Prepare the diagonal scale and validate the completed quantities
    // ---------------------------------------------------------------------
    // Guard diagonal access by the common square dimensions. A nonpositive
    // mass sum yields NaN directly, before any invalid division can occur.
    const Precision invalid    = std::numeric_limits<Precision>::quiet_NaN();
    const bool      compatible = A.rows() > 0 && A.rows() == A.cols()
                              && B.rows() == A.rows() && B.cols() == A.cols();
    const Precision stiffness_sum = compatible ? A.diagonal().sum() : invalid;
    const Precision mass_sum      = compatible ? B.diagonal().sum() : invalid;
    const Precision scale         = mass_sum > Precision(0) ? stiffness_sum / mass_sum : invalid;

    if (!compatible) return invalid;
    if (!std::isfinite(stiffness_sum) || stiffness_sum <= Precision(0)) return invalid;
    if (!std::isfinite(mass_sum) || mass_sum <= Precision(0)) return invalid;
    if (!std::isfinite(scale) || scale <= Precision(0)) return invalid;

    // ---------------------------------------------------------------------
    // (2) Return the scale in physical eigenvalue units
    // ---------------------------------------------------------------------
    // No shift factor or frequency conversion belongs to this matrix utility.
    return scale;
}

/**
 * @brief Estimates the smallest positive eigenvalue represented by quadratic fields.
 *
 * The construction uses global Cartesian nodal coordinates and the reduced
 * independent coordinates of ConstraintMap. Each of the six scalar monomials
 * x², y², z², xy, xz and yz acts in three translation directions, yielding at
 * most 18 candidate vectors. Coordinates are centered and normalized by the
 * bounding-box half-spans; derivatives are transformed back to physical units.
 * Beam/shell nodal rotations follow theta = 1/2 curl(u).
 *
 * Sampling master values defines q directly. Dependent slave values belong to
 * T q, so every trial vector satisfies the homogeneous constraints, even when
 * the original unconstrained polynomial would violate a support or coupling.
 * The affine particular solution u_p is intentionally absent from this
 * perturbation space. The procedure is not an L2 projection of a full field.
 *
 * Six global rigid-motion candidates are accepted for deflation only when
 * ||A q|| / ||A|| is small after Euclidean normalization. This is a geometric
 * candidate check, not a complete nullspace computation for mechanisms or
 * independently moving disconnected bodies. The accepted rigid basis supplies
 * B-dual vectors for deflation, with a Euclidean fallback for near-zero rigid
 * B-metric eigenvalues.
 *
 * Two-pass modified Gram-Schmidt builds a flexible basis Q with Q^T A Q close
 * to I. The projected problem is solved in reciprocal form:
 *
 *     Ar = Q^T A Q,  Br = Q^T B Q,
 *     Br a = mu Ar a,  lambda = 1 / mu.
 *
 * Ar is required to be positive definite on the accepted flexible space. Br
 * may be indefinite, so the largest positive mu selects the smallest positive
 * Ritz lambda for both modal and buckling analysis. The complete projected Ar
 * is retained to account for finite-precision orthogonalization.
 *
 * A valid FE model need not yield a useful quadratic space: all master samples
 * can vanish, a point oscillator can have zero coordinate span, or candidate
 * fields can be rigid/dependent. Such cases return NaN so the caller can choose
 * a fallback. Even a finite result may be a coarse estimate when local modes
 * are poorly represented by the global quadratic fields.
 *
 * @param A Reduced symmetric elastic stiffness, positive on flexible directions.
 * @param B Reduced symmetric mass or signed geometric-stiffness operator.
 * @param positions Global Cartesian nodal positions; only active nodes define
 *                  the bounding box, and unused node rows are ignored.
 * @param active_dof_idx_mat Node-by-six active numbering: ux, uy, uz, rx, ry, rz;
 *                           negative entries denote inactive nodal DOFs.
 * @param map Homogeneous constraint transformation with q ordered as masters.
 * @return Smallest finite positive Ritz value, or NaN for incompatible inputs,
 *         invalid arithmetic, no usable flexible basis or no positive Ritz value.
 *         All inputs and model/material states remain unchanged.
 */
Precision estimate_eigenvalue_ritz(const SparseMatrix&              A,
                                 const SparseMatrix&              B,
                                 const model::Field&              positions,
                                 const IndexMatrix&               active_dof_idx_mat,
                                 const constraint::ConstraintMap& map) {
    // ---------------------------------------------------------------------
    // (1) Validate the reduced operators and their coordinate descriptions
    // ---------------------------------------------------------------------
    // A and B must act on the same independent coordinates. The map connects
    // these coordinates to the active node/DOF numbering used below. Prepare
    // dimension-independent diagnostics first, then keep the guards together.
    const Precision invalid        = std::numeric_limits<Precision>::quiet_NaN();
    const Precision tolerance      = Precision(100) * std::numeric_limits<Precision>::epsilon();
    const Precision stiffness_norm = A.norm();

    if (A.rows() == 0 || A.rows() != A.cols() || B.rows() != A.rows() || B.cols() != A.cols()
        || map.n_master() != static_cast<Index>(A.rows()) || map.full_size != static_cast<Index>(map.T.rows())
        || map.T.cols() != A.rows() || positions.rows != static_cast<Index>(active_dof_idx_mat.rows())
        || positions.components < 3 || active_dof_idx_mat.cols() != 6) {
        return invalid;
    }
    if (!std::isfinite(stiffness_norm) || stiffness_norm <= Precision(0)) return invalid;

    // ---------------------------------------------------------------------
    // (2) Recover physical DOF identities and the active geometric extent
    // ---------------------------------------------------------------------
    // q[master] refers to the full active id map.masters[master]. Build the
    // inverse active numbering so every such id identifies a node and one of
    // its six Cartesian components. The sentinel values also detect a master
    // with no physical nodal component. Check ids before indexing the arrays.
    //
    // Only nodes with an active DOF contribute to the bounding box. Unused
    // node rows may contain NaN positions and must not distort the geometry.
    std::vector<Index> node_of_gid(map.full_size, positions.rows);
    std::vector<int>   dof_of_gid(map.full_size, -1);
    Vec3 minimum = Vec3::Constant(std::numeric_limits<Precision>::infinity());
    Vec3 maximum = -minimum;

    for (Index node = 0; node < positions.rows; ++node) {
        bool active = false;
        for (int dof = 0; dof < 6; ++dof) {
            const auto gid = active_dof_idx_mat(node, dof);
            if (gid < 0) continue; // Inactive nodal DOF.
            if (static_cast<Index>(gid) >= map.full_size || dof_of_gid[gid] >= 0) return invalid;
            node_of_gid[gid] = node;
            dof_of_gid[gid]  = dof;
            active          = true;
        }
        // A node without active DOFs has no role in this trial space.
        if (!active) continue;
        const Vec3 position(positions(node, 0), positions(node, 1), positions(node, 2));
        if (!position.allFinite()) return invalid;
        minimum = minimum.cwiseMin(position);
        maximum = maximum.cwiseMax(position);
    }
    // All subsequent trial-field loops sample only these independent DOFs.
    for (Index gid : map.masters) {
        if (gid >= map.full_size || dof_of_gid[gid] < 0) return invalid;
    }

    // ---------------------------------------------------------------------
    // (3) Define centered normalized coordinates and physical derivative scales
    // ---------------------------------------------------------------------
    // For physical X, define r = X - center and p_i = r_i / half_span_i.
    // Monomials are evaluated at dimensionless p, while gradients use
    // d/dX_i = (1 / half_span_i) d/dp_i. Physical r is retained for rotations.
    //
    // A zero span in one direction is valid for planar/linear geometries:
    // set that normalized coordinate and its derivative to zero. If all spans
    // vanish, the quadratic geometric construction is unavailable; the caller
    // can still obtain a modal shift from the separate diagonal fallback.
    const Vec3      center    = Precision(0.5) * minimum + Precision(0.5) * maximum;
    const Vec3      half_span = Precision(0.5) * maximum - Precision(0.5) * minimum;
    const Precision length    = Precision(2) * half_span.maxCoeff();
    if (!center.allFinite() || !half_span.allFinite() || !std::isfinite(length) || length <= Precision(0)) {
        return invalid;
    }
    Vec3 inverse_span = Vec3::Zero();
    for (int axis = 0; axis < 3; ++axis) {
        if (half_span[axis] > Precision(0)) inverse_span[axis] = Precision(1) / half_span[axis];
    }

    // ---------------------------------------------------------------------
    // (4) Identify independent admissible global rigid motions
    // ---------------------------------------------------------------------
    // Translation along e_i has u = e_i and theta = 0. Rotation about e_i has
    // u = e_i cross r and theta = e_i, including beam/shell rotational DOFs.
    // Restricting these values to masters can turn a physical rigid field into
    // a flexible constrained field. Therefore only a near-null vector of A is
    // eligible for removal; supported flexible content must remain available.
    std::vector<DynamicVector> rigid_basis;
    rigid_basis.reserve(6);

    for (int mode = 0; mode < 6; ++mode) {
        DynamicVector q = DynamicVector::Zero(A.rows());
        for (Index master = 0; master < map.n_master(); ++master) {
            const Index gid      = map.masters[master];
            const Index node     = node_of_gid[gid];
            const int   dof      = dof_of_gid[gid];
            const Vec3  relative(positions(node, 0) - center[0],
                                 positions(node, 1) - center[1],
                                 positions(node, 2) - center[2]);
            if (mode < 3) {
                q[master] = dof == mode ? Precision(1) : Precision(0);
            } else {
                const int axis = mode - 3;
                if (dof < 3) q[master] = Vec3::Unit(axis).cross(relative)[dof];
                else         q[master] = dof - 3 == axis ? Precision(1) : Precision(0);
            }
        }

        // Normalize before the residual test: ||A q|| / ||A|| then measures
        // stiffness response independently of the candidate's amplitude.
        const Precision norm = q.stableNorm();
        if (!std::isfinite(norm)) return invalid;
        if (norm == Precision(0)) continue; // This rigid direction has no active masters.
        q /= norm;
        const DynamicVector Aq = A * q;
        if (!Aq.allFinite()) return invalid;
        if (Aq.norm() > Precision(1e-8) * stiffness_norm) continue; // Constrained flexible field.

        // Build R with R^T R approximately I by modified Gram-Schmidt.
        // The second pass reduces loss of orthogonality for nearly dependent
        // candidates. The B metric is introduced only after R is complete;
        // using it here would assume definiteness that buckling B may not have.
        for (int pass = 0; pass < 2; ++pass) {
            for (const auto& rigid : rigid_basis) q -= rigid.dot(q) * rigid;
        }
        const Precision remaining_norm = q.stableNorm();
        if (!std::isfinite(remaining_norm)) return invalid;
        if (remaining_norm <= tolerance) continue; // Redundant rigid direction.
        q /= remaining_norm;
        rigid_basis.push_back(std::move(q));
    }

    // ---------------------------------------------------------------------
    // (5) Construct dual vectors for rigid-content removal in the B metric
    // ---------------------------------------------------------------------
    // For a nonzero modal eigenvalue and a rigid vector r with A r = 0,
    // symmetry gives 0 = r^T A q = lambda r^T B q. Flexible modes must therefore
    // be B-orthogonal to rigid motions; Euclidean orthogonality is insufficient.
    //
    // Diagonalize G = R^T B R = U diag(d_j) U^T and rotate the rigid basis:
    //     r_j = R U_j,    B r_j = B R U_j.
    // For non-negligible d_j, define dual_j = B r_j / d_j, so subtraction
    //     q <- q - r_j (dual_j^T q)
    // removes that rigid component. Both signs of d_j are admitted because B
    // may be indefinite. For near-zero d_j the existing procedure instead uses
    // dual_j = r_j, i.e. Euclidean removal of that rigid direction. In the
    // positive semidefinite modal case this includes massless rigid directions;
    // a zero projected metric alone does not prove B r_j = 0 for indefinite B.
    std::vector<DynamicVector> rigid_duals;
    if (!rigid_basis.empty()) {
        const Eigen::Index count = static_cast<Eigen::Index>(rigid_basis.size());
        DynamicMatrix      R(A.rows(), count);
        for (Eigen::Index j = 0; j < count; ++j) R.col(j) = rigid_basis[j];

        // Average with the transpose to remove floating-point asymmetry before
        // the dense symmetric eigensolve. This block has at most six rows.
        const DynamicMatrix BR           = B * R;
        DynamicMatrix       rigid_metric = R.transpose() * BR;
        rigid_metric = (Precision(0.5) * (rigid_metric + rigid_metric.transpose())).eval();
        if (!rigid_metric.allFinite()) return invalid;
        Eigen::SelfAdjointEigenSolver<DynamicMatrix> rigid_eig(rigid_metric);
        if (rigid_eig.info() != Eigen::Success) return invalid;
        const DynamicMatrix rotated_R    = R  * rigid_eig.eigenvectors();
        const DynamicMatrix rotated_BR   = BR * rigid_eig.eigenvectors();
        const Precision     metric_scale = rigid_eig.eigenvalues().cwiseAbs().maxCoeff();
        for (Eigen::Index j = 0; j < count; ++j) {
            rigid_basis[j] = rotated_R.col(j);
            const Precision value = rigid_eig.eigenvalues()[j];
            if (std::abs(value) > tolerance * metric_scale) {
                rigid_duals.push_back(rotated_BR.col(j) / value);
            } else {
                rigid_duals.push_back(rotated_R.col(j));
            }
        }
    }

    // ---------------------------------------------------------------------
    // (6) Sample quadratic fields and build the flexible elastic-energy basis
    // ---------------------------------------------------------------------
    // Six monomials times three translation directions produce at most 18
    // candidates. For monomial f(p), direction e_i and model length L, use
    //     u(X)     = L f(p(X)) e_i,
    //     theta(X) = (L / 2) grad_X(f) cross e_i.
    // L is the maximum bounding-box span. It gives translations length units
    // and rotations angle units, preserving their compatible relative scale.
    //
    // Values are assigned only to masters. Slave values are implicitly T q;
    // they need not equal the unconstrained polynomial at their own nodes.
    // stiffness_basis[j] stores A * ritz_basis[j] to reuse the sparse product
    // during elastic-energy Gram-Schmidt.
    std::vector<DynamicVector> ritz_basis;
    std::vector<DynamicVector> stiffness_basis;
    ritz_basis.reserve(18);
    stiffness_basis.reserve(18);

    for (int polynomial = 0; polynomial < 6; ++polynomial) {
        for (int direction = 0; direction < 3; ++direction) {
            DynamicVector q = DynamicVector::Zero(A.rows());
            for (Index master = 0; master < map.n_master(); ++master) {
                const Index gid      = map.masters[master];
                const Index node     = node_of_gid[gid];
                const int   dof      = dof_of_gid[gid];
                const Vec3  relative(positions(node, 0) - center[0],
                                     positions(node, 1) - center[1],
                                     positions(node, 2) - center[2]);
                const Vec3  p        = relative.cwiseProduct(inverse_span);

                // Scalar fields use normalized coordinates. Their gradients
                // include inverse_span so curl(u) uses physical derivatives.
                const std::array<Precision, 6> monomials{
                    p[0] * p[0], p[1] * p[1], p[2] * p[2],
                    p[0] * p[1], p[0] * p[2], p[1] * p[2]
                };
                const std::array<Vec3, 6> gradients{
                    Vec3(Precision(2) * p[0] * inverse_span[0], 0, 0),
                    Vec3(0, Precision(2) * p[1] * inverse_span[1], 0),
                    Vec3(0, 0, Precision(2) * p[2] * inverse_span[2]),
                    Vec3(p[1] * inverse_span[0], p[0] * inverse_span[1], 0),
                    Vec3(p[2] * inverse_span[0], 0,                       p[0] * inverse_span[2]),
                    Vec3(0,                       p[2] * inverse_span[1], p[1] * inverse_span[2])
                };

                // Translational and rotational samples describe the same
                // physical displacement field, rather than independent probes.
                if (dof == direction) q[master] = length * monomials[polynomial];
                if (dof >= 3) {
                    const Vec3 rotation = Precision(0.5) * length * gradients[polynomial].cross(Vec3::Unit(direction));
                    q[master] = rotation[dof - 3];
                }
            }

            // First remove amplitude dependence, then remove rigid content in
            // two passes. A geometric zero sample is a valid empty candidate,
            // for example when all translational masters lie at the center.
            const Precision initial_norm = q.stableNorm();
            if (!q.allFinite() || !std::isfinite(initial_norm)) return invalid;
            if (initial_norm == Precision(0)) continue; // Vanishing monomial on this geometry.
            q /= initial_norm;
            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t j = 0; j < rigid_basis.size(); ++j) {
                    q -= rigid_duals[j].dot(q) * rigid_basis[j];
                }
            }

            // Renormalize after deflation so the subsequent stiffness-energy
            // threshold is applied to an initially Euclidean unit direction.
            const Precision deflated_norm = q.stableNorm();
            if (!q.allFinite() || !std::isfinite(deflated_norm)) return invalid;
            if (deflated_norm <= tolerance) continue; // Pure rigid trial field.
            q /= deflated_norm;

            // Modified Gram-Schmidt in <x,y>_A = x^T A y:
            //     alpha = Q_j^T A q,
            //     q    <- q    - alpha Q_j,
            //     A q  <- A q  - alpha A Q_j.
            // Accepted Q_j have unit elastic energy. Updating q and A q in
            // parallel avoids another sparse matrix-vector product per step.
            // Two passes reduce cancellation errors for dependent monomials;
            // using A avoids a positive-metric assumption on buckling B.
            DynamicVector Aq = A * q;
            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t j = 0; j < ritz_basis.size(); ++j) {
                    const Precision projection = ritz_basis[j].dot(Aq);
                    q  -= projection * ritz_basis[j];
                    Aq -= projection * stiffness_basis[j];
                }
            }
            const Precision energy = q.dot(Aq);
            if (!q.allFinite() || !Aq.allFinite()) return invalid;
            if (!std::isfinite(energy) || energy < -tolerance * stiffness_norm) return invalid;

            // A significantly negative energy violates the assumed flexible
            // stiffness metric. A zero/tiny nonnegative energy instead marks a
            // null or nearly dependent candidate and is skipped. The threshold
            // uses ||A|| and the unit amplitude before Gram-Schmidt, rather than
            // the norm of the small remainder left after cancellation.
            if (energy <= tolerance * stiffness_norm) continue;

            // Scale q so q^T A q = 1, retaining its corresponding A q.
            const Precision inverse_norm = Precision(1) / std::sqrt(energy);
            q  *= inverse_norm;
            Aq *= inverse_norm;
            ritz_basis.push_back(std::move(q));
            stiffness_basis.push_back(std::move(Aq));
        }
    }
    // No flexible field can be inferred from an empty space. This says nothing
    // about whether the full FE model has positive eigenvalues; the modal
    // caller can still use the independent diagonal-sum fallback.
    if (ritz_basis.empty()) return invalid;

    // ---------------------------------------------------------------------
    // (7) Project both operators onto the accepted quadratic subspace
    // ---------------------------------------------------------------------
    // For a reduced trial vector phi = Q a, Galerkin projection gives
    //     Ar a = lambda Br a,   Ar = Q^T A Q,   Br = Q^T B Q.
    // Recompute A Q from the final Q to make the projected operator consistent
    // despite roundoff accumulated in the cached Gram-Schmidt products. Ar is
    // close to I, but explicitly forming it preserves the actual Ritz problem.
    const Eigen::Index count = static_cast<Eigen::Index>(ritz_basis.size());
    DynamicMatrix      Q(A.rows(), count);
    for (Eigen::Index j = 0; j < count; ++j) {
        Q.col(j) = ritz_basis[j];
    }
    const DynamicMatrix AQ          = A * Q;
    const DynamicMatrix BQ          = B * Q;
    DynamicMatrix       projected_A = Q.transpose() * AQ;
    DynamicMatrix       projected_B = Q.transpose() * BQ;

    // Both reduced operators are symmetric by construction. Average paired
    // entries before passing the small dense matrices to the symmetric solver.
    projected_A = (Precision(0.5) * (projected_A + projected_A.transpose())).eval();
    projected_B = (Precision(0.5) * (projected_B + projected_B.transpose())).eval();
    if (!projected_A.allFinite() || !projected_B.allFinite()) return invalid;

    // ---------------------------------------------------------------------
    // (8) Solve the reciprocal problem and select the first positive Ritz value
    // ---------------------------------------------------------------------
    // The generalized self-adjoint solver needs a positive definite right-hand
    // metric. Solve Br a = mu Ar a so Ar supplies that metric even for an
    // indefinite Br. For positive values, lambda = 1 / mu reverses the order:
    // the largest positive mu corresponds to the smallest positive lambda.
    // Eigen returns the reciprocal eigenvalues in ascending order.
    Eigen::GeneralizedSelfAdjointEigenSolver<DynamicMatrix> eig(projected_B, projected_A, Eigen::EigenvaluesOnly);
    if (eig.info() != Eigen::Success) return invalid;
    const Precision reciprocal = eig.eigenvalues()[count - 1];
    if (!std::isfinite(reciprocal) || reciprocal <= Precision(0)) return invalid;
    // Keep the result in the caller's original eigenvalue units. A failed or
    // nonpositive reciprocal cannot define a usable spectral scale here.
    const Precision estimate = Precision(1) / reciprocal;
    return std::isfinite(estimate) && estimate > Precision(0) ? estimate : invalid;
}

} // namespace mattools
} // namespace fem
