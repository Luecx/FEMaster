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

#include "../core/logging.h"

#include <array>
#include <cmath>
#include <limits>
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
 * elastic stiffness and mass. Zero stiffness gives a valid zero scale; the
 * caller supplies a nonzero shift reference for purely zero-frequency modes.
 * The sums are formed before division, so a DOF
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
 * @return Finite nonnegative diagonal-sum ratio, or NaN for incompatible dimensions,
 *         nonfinite arithmetic, negative stiffness sum or nonpositive mass sum. No input is
 *         modified and no factorization is performed.
 */
Precision estimate_spectral_scale(const SparseMatrix& A, const SparseMatrix& B) {
    // ---------------------------------------------------------------------
    // (1) Prepare the diagonal scale and validate the completed quantities
    // ---------------------------------------------------------------------
    // First verify matching square dimensions, then sum the diagonals.
    // Forming both sums before division avoids division by a zero mass entry.
    const Precision invalid = std::numeric_limits<Precision>::quiet_NaN();
    if (A.rows() <= 0 || A.rows() != A.cols()) return invalid;
    if (B.rows() != A.rows() || B.cols() != A.cols()) return invalid;

    const Precision stiffness_sum = A.diagonal().sum();
    const Precision mass_sum      = B.diagonal().sum();
    if (!std::isfinite(stiffness_sum) || stiffness_sum < 0) return invalid;
    if (!std::isfinite(mass_sum) || mass_sum <= 0) return invalid;

    const Precision scale = stiffness_sum / mass_sum;
    if (!std::isfinite(scale) || scale < 0) return invalid;

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
Precision estimate_eigenvalue_ritz(
     const SparseMatrix&              A,
     const SparseMatrix&              B,
     const model::Field&              positions,
     const IndexMatrix&               active_dof_idx_mat,
     const constraint::ConstraintMap& map) {
    // ---------------------------------------------------------------------
    // (1) Validate the reduced operators and their coordinate descriptions
    // ---------------------------------------------------------------------
    // A and B act on the same independent coordinates. The map connects these
    // coordinates to the full active node/DOF numbering. Structural mismatches
    // are programming errors; an empty or numerically unusable operator may
    // legitimately fall back to a different spectral scale estimator.
    const Precision invalid        = std::numeric_limits<Precision>::quiet_NaN();
    const Precision tolerance      = Precision(100) * std::numeric_limits<Precision>::epsilon();
    const Precision stiffness_norm = A.norm();

    // Validate matrix dimensions and constraint mapping
    logging::error(A.rows()  == A.cols(),
        "Ritz estimator requires a square stiffness matrix");
    logging::error(B.rows()  == A.rows() && B.cols() == A.cols(),
        "Ritz estimator requires matching operator dimensions");
    logging::error(map.n_master() == static_cast<Index>(A.rows()),
        "Ritz estimator: master count does not match reduced system");
    logging::error(map.full_size == static_cast<Index>(map.T.rows()),
        "Ritz estimator: invalid constraint transformation rows");
    logging::error(map.T.cols() == A.rows(),
        "Ritz estimator: invalid constraint transformation columns");

    // Validate nodal geometry and DOF layout
    logging::error(positions.rows == static_cast<Index>(active_dof_idx_mat.rows()),
        "Ritz estimator: nodal position and DOF counts differ");
    logging::error(positions.components >= 3,
        "Ritz estimator requires three-dimensional positions");
    logging::error(active_dof_idx_mat.cols() == 6,
        "Ritz estimator requires six DOF columns");

    // Reject systems without a usable elastic stiffness scale
    if (A.rows() == 0) return invalid;
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
    // Compute the bounding-box center, half-spans and characteristic length
    const Vec3      center    = minimum / 2 + maximum / 2;
    const Vec3      half_span = maximum / 2 - minimum / 2;
    const Precision length    = 2 * half_span.maxCoeff();

    if (!std::isfinite(length) || length <= 0) return invalid;

    // Compute inverse half-spans, preserving zero for degenerate dimensions
    Vec3 inverse_span = Vec3::Zero();
    for (int axis = 0; axis < 3; ++axis) {
        if (half_span[axis] > 0)
            inverse_span[axis] = Precision(1) / half_span[axis];
    }

    // ---------------------------------------------------------------------
    // (4) Identify independent admissible global rigid motions
    // ---------------------------------------------------------------------
    // For each Cartesian axis e_i, a global translation has u = e_i and
    // theta = 0; a global rotation has u = e_i cross (X - center) and
    // theta = e_i. The first three columns of R represent translations,
    // followed by the three rotations. Sampling only masters gives reduced
    // candidates; their complete constrained fields would be T * R.
    //
    // Such sampled motions need not be truly rigid after constraint recovery:
    // an MPC can map their master values to non-rigid slave displacements.
    // Therefore accept a candidate for deflation only if A q is near zero.
    // R is reused for the accepted orthonormal basis, avoiding a second
    // collection of full-length vectors.
    DynamicMatrix R = DynamicMatrix::Zero(A.rows(), 6);

    for (Index master = 0; master < map.n_master(); ++master) {
        const Index gid  = map.masters[master];
        const int   dof  = dof_of_gid[gid];

        if (dof < 3) {
            const Index node = node_of_gid[gid];
            const Vec3  r(positions(node, 0) - center[0],
                          positions(node, 1) - center[1],
                          positions(node, 2) - center[2]);

            R(master, dof) = 1;
            for (int axis = 0; axis < 3; ++axis) {
                R(master, axis + 3) = Vec3::Unit(axis).cross(r)[dof];
            }
        } else {
            // Rotational DOFs have no translational rigid-mode component.
            R(master, dof) = 1;
        }
    }

    // Test the original six candidates in column order. A rigid displacement
    // must satisfy A q approximately zero; normalize first so the residual
    // threshold is independent of the candidate's arbitrary amplitude.
    // Use two-pass Euclidean modified Gram-Schmidt to remove dependent rigid
    // candidates. The B metric is not suitable here because it may be
    // indefinite in buckling analyses.
    Eigen::Index rigid_count = 0;
    for (int mode = 0; mode < 6; ++mode) {
        DynamicVector q = R.col(mode);
        const Precision norm = q.stableNorm();
        if (!std::isfinite(norm)) return invalid;
        if (norm == 0) continue; // This mode has no active master component.
        q /= norm;

        const DynamicVector Aq = A * q;
        if (!Aq.allFinite()) return invalid;
        if (Aq.norm() > Precision(1e-8) * stiffness_norm) continue; // Not a free rigid mode.

        for (int pass = 0; pass < 2; ++pass) {
            for (Eigen::Index j = 0; j < rigid_count; ++j) {
                q -= R.col(j).dot(q) * R.col(j);
            }
        }

        const Precision remaining_norm = q.stableNorm();
        if (!std::isfinite(remaining_norm)) return invalid;
        if (remaining_norm <= tolerance) continue; // Linearly dependent candidate.
        R.col(rigid_count++) = q / remaining_norm;
    }
    R.conservativeResize(Eigen::NoChange, rigid_count);

    // ---------------------------------------------------------------------
    // (5) Construct B-dual vectors for rigid-content removal
    // ---------------------------------------------------------------------
    // For a nonzero eigenvalue and a rigid vector r with A r = 0, symmetry
    // gives r^T A q = lambda r^T B q = 0. Thus flexible eigenvectors are
    // B-orthogonal to the stiffness nullspace, not generally Euclidean-
    // orthogonal. We remove the B component along the accepted rigid modes.
    //
    // Diagonalize the small symmetric rigid metric
    //     G = R^T B R = U diag(d_j) U^T.
    // Rotate the rigid columns to R U and form corresponding dual columns
    //     D_j = B (R U)_j / d_j   for sufficiently large |d_j|.
    // This supports positive and negative metric eigenvalues. For near-zero
    // d_j, use D_j = (R U)_j instead, producing Euclidean deflation in that
    // direction. With indefinite B, a zero d_j does not imply B R_j = 0.
    // In both cases D^T R is approximately the identity, so the two-pass
    // subtraction q <- q - R (D^T q) removes rigid content.
    DynamicMatrix rigid_duals(A.rows(), rigid_count);
    if (rigid_count > 0) {
        const DynamicMatrix BR = B * R;
        DynamicMatrix rigid_metric = R.transpose() * BR;
        rigid_metric = (Precision(0.5) * (rigid_metric + rigid_metric.transpose())).eval();
        if (!rigid_metric.allFinite()) return invalid;

        Eigen::SelfAdjointEigenSolver<DynamicMatrix> rigid_eig(rigid_metric);
        if (rigid_eig.info() != Eigen::Success) return invalid;

        // Rotate both R and B R with the same eigenvectors. Only six or fewer
        // columns are present; this does not factorize the full FE operators.
        const DynamicMatrix rotated_R = R * rigid_eig.eigenvectors();
        rigid_duals = BR * rigid_eig.eigenvectors();
        R = rotated_R;

        const Precision metric_scale = rigid_eig.eigenvalues().cwiseAbs().maxCoeff();
        for (Eigen::Index j = 0; j < rigid_count; ++j) {
            const Precision value = rigid_eig.eigenvalues()[j];
            if (std::abs(value) > tolerance * metric_scale) {
                rigid_duals.col(j) /= value;
            } else {
                rigid_duals.col(j) = R.col(j);
            }
        }
    }

    // ---------------------------------------------------------------------
    // (6) Sample quadratic fields and build the flexible elastic-energy basis
    // ---------------------------------------------------------------------
    // Each of the six scalar monomials f(p) = x², y², z², xy, xz, yz is
    // combined with all three Cartesian displacement directions e_i:
    //     u(X)     = L f(p(X)) e_i,
    //     theta(X) = (L / 2) grad_X(f) cross e_i.
    // The characteristic length L scales translation and rotation consistently.
    // Derivatives contain inverse half-spans because the monomials use
    // dimensionless coordinates. Slave values are not sampled separately;
    // they are implied by T q and need not match the unconstrained polynomial.
    //
    // For each monomial, construct the three displacement directions together.
    // This evaluates coordinates, monomials and gradients once per master,
    // rather than repeating the same computation for each direction. Only
    // three temporary candidate columns are needed at any time.
    //
    // Q stores the accepted elastic-energy-orthonormal Ritz vectors, while KQ
    // caches A * Q for two-pass modified Gram-Schmidt. Both have up to 18
    // columns; the current count identifies the accepted leading columns.
    DynamicMatrix Q(A.rows(), 18);
    DynamicMatrix KQ(A.rows(), 18);
    Eigen::Index  count = 0;

    for (int polynomial = 0; polynomial < 6; ++polynomial) {
        DynamicMatrix candidates = DynamicMatrix::Zero(A.rows(), 3);

        for (Index master = 0; master < map.n_master(); ++master) {
            const Index gid      = map.masters[master];
            const Index node     = node_of_gid[gid];
            const int   dof      = dof_of_gid[gid];
            const Vec3  relative(positions(node, 0) - center[0],
                                 positions(node, 1) - center[1],
                                 positions(node, 2) - center[2]);
            const Vec3 p = relative.cwiseProduct(inverse_span);

            // Scalar polynomial values are dimensionless. Gradients are with
            // respect to physical X, not p, for the compatible nodal curl.
            const Precision x  = p[0], y  = p[1], z  = p[2];
            const Precision ix = inverse_span[0], iy = inverse_span[1], iz = inverse_span[2];

            const std::array<Precision, 6> monomials{
                x * x, y * y, z * z, x * y, x * z, y * z
            };
            const std::array<Vec3, 6> gradients{
                Vec3(2 * x * ix, 0,          0),
                Vec3(0,          2 * y * iy, 0),
                Vec3(0,          0,          2 * z * iz),
                Vec3(y * ix,     x * iy,     0),
                Vec3(z * ix,     0,          x * iz),
                Vec3(0,          z * iy,     y * iz)
            };

            if (dof < 3) {
                // Only displacement in the selected direction is nonzero.
                candidates(master, dof) = length * monomials[polynomial];
            } else {
                // Nodal rotation equals half the curl of the same field.
                for (int direction = 0; direction < 3; ++direction) {
                    candidates(master, direction) = Precision(0.5) * length
                        * gradients[polynomial].cross(Vec3::Unit(direction))[dof - 3];
                }
            }
        }

        // Treat the three fields in the same order as the original algorithm:
        // first polynomial, then x/y/z direction. Preserve this order for
        // reproducible orthogonalization and dependence checks.
        for (int direction = 0; direction < 3; ++direction) {
            DynamicVector q = candidates.col(direction);

            // Normalize first and remove rigid content in two passes. A zero
            // sample is valid (e.g. a monomial vanishing on all masters).
            const Precision initial_norm = q.stableNorm();
            if (!q.allFinite() || !std::isfinite(initial_norm)) return invalid;
            if (initial_norm == 0) continue;
            q /= initial_norm;

            for (int pass = 0; pass < 2; ++pass) {
                if (rigid_count == 0) break;
                const DynamicVector weights = rigid_duals.transpose() * q;
                q -= R * weights;
            }

            // Normalize after deflation so the stiffness-energy threshold is
            // relative to a Euclidean unit vector, not the arbitrary amplitude
            // left after subtracting rigid content.
            const Precision deflated_norm = q.stableNorm();
            if (!q.allFinite() || !std::isfinite(deflated_norm)) return invalid;
            if (deflated_norm <= tolerance) continue; // Pure rigid direction.
            q /= deflated_norm;

            // Two-pass modified Gram-Schmidt in the elastic energy metric:
            //     alpha = Q_j^T A q,
            //     q     <- q   - alpha Q_j,
            //     A q   <- A q - alpha A Q_j.
            // Reuse KQ rather than recomputing a sparse product for every
            // projection. B may be indefinite, so it is not a valid positive
            // orthogonalization metric in the buckling case.
            DynamicVector Aq = A * q;
            for (int pass = 0; pass < 2; ++pass) {
                for (Eigen::Index j = 0; j < count; ++j) {
                    const Precision projection = Q.col(j).dot(Aq);
                    q  -= projection * Q.col(j);
                    Aq -= projection * KQ.col(j);
                }
            }

            // A is assumed nonnegative on flexible directions. Reject negative
            // energy beyond roundoff; zero/tiny energy means a null or dependent
            // candidate. Compare against ||A|| at unit preprojection amplitude.
            const Precision energy = q.dot(Aq);
            if (!q.allFinite() || !Aq.allFinite()) return invalid;
            if (!std::isfinite(energy) || energy < -tolerance * stiffness_norm) return invalid;
            if (energy <= tolerance * stiffness_norm) continue;

            // Store Q_j^T A Q_j = 1 together with the corresponding A Q_j.
            const Precision inverse_norm = Precision(1) / std::sqrt(energy);
            Q.col(count)  = inverse_norm * q;
            KQ.col(count) = inverse_norm * Aq;
            ++count;
        }
    }

    // An empty quadratic flexible space does not imply that the FE model has
    // no positive eigenvalues. The caller may use its independent fallback.
    if (count == 0) return invalid;
    Q.conservativeResize(Eigen::NoChange, count);

    // ---------------------------------------------------------------------
    // (7) Project both operators onto the accepted quadratic subspace
    // ---------------------------------------------------------------------
    // For a reduced trial vector phi = Q a, Galerkin projection gives
    //     Ar a = lambda Br a,   Ar = Q^T A Q,   Br = Q^T B Q.
    // Recompute A Q from the accepted Q, rather than relying on KQ after
    // repeated orthogonalization updates. Ar is close to I, but retaining it
    // explicitly preserves the actual projected Ritz problem in finite precision.
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
