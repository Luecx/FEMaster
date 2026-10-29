/**
 * @file eigenvalue_estimates.cpp
 * @brief Implements quadratic Rayleigh-Ritz estimates for structural analyses.
 *
 * Normalized nodal coordinates generate quadratic displacement fields and
 * compatible rotations. Independent master values define admissible trial
 * vectors through ConstraintMap; admissible rigid motions are removed before
 * elastic-energy orthonormalization. A small symmetric reciprocal eigenproblem
 * handles both positive mass and indefinite geometric stiffness matrices.
 * A diagonal-sum scale supplies a geometry-independent modal fallback when
 * the quadratic trial space contains no usable flexible direction.
 *
 * The calling load case assembles the operators, supplies the sign of geometric
 * stiffness and chooses its spectral shift. No full-system factorization is
 * performed and no input or material state is changed.
 *
 * @see eigenvalue_estimates.h
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
 * @brief Returns a characteristic modal eigenvalue scale from diagonal sums.
 *
 * For positive semidefinite stiffness A and mass B, the ratio trace(A)/trace(B)
 * has the units of omega^2 and is positive whenever elastic stiffness and mass
 * are present. Summing before division also supports massless individual DOFs.
 * Geometry, coordinate spans and availability of polynomial fields play no role.
 * This scale does not bound or specifically estimate the lowest eigenvalue.
 *
 * The modal caller selects this fallback only when its Ritz estimate is invalid
 * and retains ownership of the shift. Indefinite geometric-stiffness operators
 * should continue to use the Ritz estimator rather than this modal scale.
 *
 * @return Finite positive diagonal-sum ratio, or NaN for incompatible matrices,
 *         nonfinite arithmetic or nonpositive stiffness/mass sums.
 */
Precision estimate_spectral_scale(const SparseMatrix& A, const SparseMatrix& B) {
    // Prepare safe diagnostics first, then validate them as one consecutive block.
    const Precision invalid = std::numeric_limits<Precision>::quiet_NaN();
    const bool compatible = A.rows() > 0 && A.rows() == A.cols()
                         && B.rows() == A.rows() && B.cols() == A.cols();
    const Precision stiffness_sum = compatible ? A.diagonal().sum() : invalid;
    const Precision mass_sum      = compatible ? B.diagonal().sum() : invalid;
    const Precision scale         = mass_sum > Precision(0) ? stiffness_sum / mass_sum : invalid;

    if (!compatible) return invalid;
    if (!std::isfinite(stiffness_sum) || stiffness_sum <= Precision(0)) return invalid;
    if (!std::isfinite(mass_sum) || mass_sum <= Precision(0)) return invalid;
    if (!std::isfinite(scale) || scale <= Precision(0)) return invalid;

    return scale;
}

/**
 * @brief Estimates the lowest positive eigenvalue in a quadratic Ritz space.
 *
 * The six monomials x², y², z², xy, xz and yz use centered coordinates scaled
 * by each geometric half-span. Each monomial acts in all three global
 * translation directions; nodal rotations follow theta = 1/2 curl(u), with
 * derivatives taken in physical coordinates. Sampling only independent master
 * DOFs makes the recovered fields T q satisfy homogeneous constraints.
 *
 * Rigid translations and rotations are retained for deflation only when they
 * lie numerically in the nullspace of the reduced stiffness A. Modified
 * Gram-Schmidt then builds an A-orthonormal flexible basis Q. The projected
 * problem is solved in reciprocal form
 *
 *     (Q^T B Q) a = mu (Q^T A Q) a,  lambda = 1 / mu.
 *
 * Q^T A Q is positive definite, so B need not be positive definite. The largest
 * positive mu gives the smallest positive Ritz lambda, for both B = M and
 * B = -Kg. Vanishing coordinate spans, redundant trial fields and rigid fields
 * do not contribute to the flexible basis. Sparse unused node rows are ignored.
 *
 * @return Smallest finite positive Ritz value, or NaN for incompatible input,
 *         invalid arithmetic, an empty flexible space or no positive Ritz value.
 */
Precision estimate_eigenvalue_ritz(const SparseMatrix&              A,
                                   const SparseMatrix&              B,
                                   const model::Field&              positions,
                                   const IndexMatrix&               active_dof_idx_mat,
                                   const constraint::ConstraintMap& map) {
    // Validate the common reduced coordinate space before indexing or products.
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

    // Map active system indices back to node/DOF pairs. Only active nodes enter
    // the bounding box; unassigned node rows may contain nonfinite coordinates.
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
            active = true;
        }
        if (!active) continue;
        const Vec3 position(positions(node, 0), positions(node, 1), positions(node, 2));
        if (!position.allFinite()) return invalid;
        minimum = minimum.cwiseMin(position);
        maximum = maximum.cwiseMax(position);
    }
    for (Index gid : map.masters) {
        if (gid >= map.full_size || dof_of_gid[gid] < 0) return invalid;
    }

    // Normalize coordinates to avoid powers of the absolute model location.
    // Zero spans describe planar/linear models and yield zero derivatives.
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

    // Detect admissible rigid motions in reduced coordinates. Restricting a
    // physical rigid motion to masters alone can alter it at supported/slave
    // DOFs; the stiffness residual therefore decides whether it is removable.
    std::vector<DynamicVector> rigid_basis;
    rigid_basis.reserve(6);
    for (int mode = 0; mode < 6; ++mode) {
        DynamicVector q = DynamicVector::Zero(A.rows());
        for (Index master = 0; master < map.n_master(); ++master) {
            const Index gid  = map.masters[master];
            const Index node = node_of_gid[gid];
            const int   dof  = dof_of_gid[gid];
            const Vec3 relative(positions(node, 0) - center[0],
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
        const Precision norm = q.stableNorm();
        if (!std::isfinite(norm)) return invalid;
        if (norm == Precision(0)) continue; // This rigid direction has no active masters.
        q /= norm;
        const DynamicVector Aq = A * q;
        if (!Aq.allFinite()) return invalid;
        if (Aq.norm() > Precision(1e-8) * stiffness_norm) continue; // Constrained flexible field.

        // Build an independent Euclidean basis first; the B metric used to
        // deflate trial fields is introduced below after all rigid modes exist.
        for (int pass = 0; pass < 2; ++pass) {
            for (const auto& rigid : rigid_basis) q -= rigid.dot(q) * rigid;
        }
        const Precision remaining_norm = q.stableNorm();
        if (!std::isfinite(remaining_norm)) return invalid;
        if (remaining_norm <= tolerance) continue; // Redundant rigid direction.
        q /= remaining_norm;
        rigid_basis.push_back(std::move(q));
    }

    // Flexible vibration modes are B-orthogonal to rigid motions, rather
    // than Euclidean-orthogonal. Diagonalize the small rigid B block to obtain
    // dual vectors for that projection without assuming B is positive definite.
    // For massless rigid directions, retain Euclidean deflation; their B norm
    // vanishes and they cannot supply a finite flexible eigenvalue.
    std::vector<DynamicVector> rigid_duals;
    if (!rigid_basis.empty()) {
        const Eigen::Index count = static_cast<Eigen::Index>(rigid_basis.size());
        DynamicMatrix R(A.rows(), count);
        for (Eigen::Index j = 0; j < count; ++j) R.col(j) = rigid_basis[j];
        const DynamicMatrix BR = B * R;
        DynamicMatrix rigid_metric = R.transpose() * BR;
        rigid_metric = (Precision(0.5) * (rigid_metric + rigid_metric.transpose())).eval();
        if (!rigid_metric.allFinite()) return invalid;
        Eigen::SelfAdjointEigenSolver<DynamicMatrix> rigid_eig(rigid_metric);
        if (rigid_eig.info() != Eigen::Success) return invalid;
        const DynamicMatrix rotated_R  = R  * rigid_eig.eigenvectors();
        const DynamicMatrix rotated_BR = BR * rigid_eig.eigenvectors();
        const Precision metric_scale = rigid_eig.eigenvalues().cwiseAbs().maxCoeff();
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

    // Generate quadratic translations and their physical curls directly at
    // masters. Recovered slaves follow T, rather than an unconstrained field.
    // The factor L keeps translation amplitudes comparable to nodal rotations.
    std::vector<DynamicVector> ritz_basis;
    std::vector<DynamicVector> stiffness_basis;
    ritz_basis.reserve(18);
    stiffness_basis.reserve(18);
    for (int polynomial = 0; polynomial < 6; ++polynomial) {
        for (int direction = 0; direction < 3; ++direction) {
            DynamicVector q = DynamicVector::Zero(A.rows());
            for (Index master = 0; master < map.n_master(); ++master) {
                const Index gid  = map.masters[master];
                const Index node = node_of_gid[gid];
                const int   dof  = dof_of_gid[gid];
                const Vec3 relative(positions(node, 0) - center[0],
                                    positions(node, 1) - center[1],
                                    positions(node, 2) - center[2]);
                const Vec3 p = relative.cwiseProduct(inverse_span);
                const std::array<Precision, 6> monomials{
                    p[0] * p[0], p[1] * p[1], p[2] * p[2],
                    p[0] * p[1], p[0] * p[2], p[1] * p[2]
                };
                const std::array<Vec3, 6> gradients{
                    Vec3(Precision(2) * p[0] * inverse_span[0], 0, 0),
                    Vec3(0, Precision(2) * p[1] * inverse_span[1], 0),
                    Vec3(0, 0, Precision(2) * p[2] * inverse_span[2]),
                    Vec3(p[1] * inverse_span[0], p[0] * inverse_span[1], 0),
                    Vec3(p[2] * inverse_span[0], 0, p[0] * inverse_span[2]),
                    Vec3(0, p[2] * inverse_span[1], p[1] * inverse_span[2])
                };
                if (dof == direction) q[master] = length * monomials[polynomial];
                if (dof >= 3) {
                    const Vec3 rotation = Precision(0.5) * length * gradients[polynomial].cross(Vec3::Unit(direction));
                    q[master] = rotation[dof - 3];
                }
            }
            const Precision initial_norm = q.stableNorm();
            if (!q.allFinite() || !std::isfinite(initial_norm)) return invalid;
            if (initial_norm == Precision(0)) continue; // Vanishing monomial on this geometry.
            q /= initial_norm;
            for (int pass = 0; pass < 2; ++pass) {
                for (std::size_t j = 0; j < rigid_basis.size(); ++j) {
                    q -= rigid_duals[j].dot(q) * rigid_basis[j];
                }
            }
            const Precision deflated_norm = q.stableNorm();
            if (!q.allFinite() || !std::isfinite(deflated_norm)) return invalid;
            if (deflated_norm <= tolerance) continue; // Pure rigid trial field.
            q /= deflated_norm;

            // Orthogonalize in elastic energy, never in the potentially
            // indefinite geometric-stiffness metric. Update A q alongside q.
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
            // Reject null/near-dependent fields using the stiffness scale of
            // the initially unit-normalized candidate, before cancellation.
            if (energy <= tolerance * stiffness_norm) continue;
            const Precision inverse_norm = Precision(1) / std::sqrt(energy);
            q  *= inverse_norm;
            Aq *= inverse_norm;
            ritz_basis.push_back(std::move(q));
            stiffness_basis.push_back(std::move(Aq));
        }
    }
    if (ritz_basis.empty()) return invalid;

    // Form the small reciprocal problem. Retaining the actual projected A
    // accounts for finite-precision orthogonalization instead of assuming I.
    const Eigen::Index count = static_cast<Eigen::Index>(ritz_basis.size());
    DynamicMatrix Q(A.rows(), count);
    for (Eigen::Index j = 0; j < count; ++j) {
        Q.col(j) = ritz_basis[j];
    }
    const DynamicMatrix AQ = A * Q;
    const DynamicMatrix BQ = B * Q;
    DynamicMatrix projected_A = Q.transpose() * AQ;
    DynamicMatrix projected_B = Q.transpose() * BQ;
    projected_A = (Precision(0.5) * (projected_A + projected_A.transpose())).eval();
    projected_B = (Precision(0.5) * (projected_B + projected_B.transpose())).eval();
    if (!projected_A.allFinite() || !projected_B.allFinite()) return invalid;

    Eigen::GeneralizedSelfAdjointEigenSolver<DynamicMatrix> eig(projected_B, projected_A, Eigen::EigenvaluesOnly);
    if (eig.info() != Eigen::Success) return invalid;
    const Precision reciprocal = eig.eigenvalues()[count - 1];
    if (!std::isfinite(reciprocal) || reciprocal <= Precision(0)) return invalid;
    const Precision estimate = Precision(1) / reciprocal;
    return std::isfinite(estimate) && estimate > Precision(0) ? estimate : invalid;
}

} // namespace mattools
} // namespace fem
