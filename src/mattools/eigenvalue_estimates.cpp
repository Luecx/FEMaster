/**
 * @file eigenvalue_estimates.cpp
 * @brief Estimates modal and buckling eigenvalues without factorizing the FE system.
 *
 * @author Finn Eggers
 * @date 29.10.2026
 */

#include "eigenvalue_estimates.h"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>
#include <Eigen/SVD>

#include "../core/logging.h"

#include <array>
#include <cmath>
#include <limits>

namespace fem {
namespace mattools {

namespace {

// Remove small asymmetries introduced by floating-point matrix products.
DynamicMatrix symmetric(const DynamicMatrix& A) {
    return (Precision(0.5) * (A + A.transpose())).eval();
}

} // namespace

/**
 * @brief Returns a characteristic modal eigenvalue scale, trace(A) / trace(B).
 * @return Nonnegative scale, or NaN if the matrices are incompatible or invalid.
 */
Precision estimate_spectral_scale(const SparseMatrix& A, const SparseMatrix& B) {
    const Precision invalid = std::numeric_limits<Precision>::quiet_NaN();
    if (A.rows() <= 0 || A.rows() != A.cols()) return invalid;
    if (B.rows() != A.rows() || B.cols() != A.cols()) return invalid;

    const Precision stiffness_sum = A.diagonal().sum();
    const Precision mass_sum      = B.diagonal().sum();
    if (!std::isfinite(stiffness_sum) || stiffness_sum < 0) return invalid;
    if (!std::isfinite(mass_sum) || mass_sum <= 0) return invalid;

    const Precision scale = stiffness_sum / mass_sum;
    return std::isfinite(scale) && scale >= 0 ? scale : invalid;
}

/**
 * @brief Estimates the smallest positive eigenvalue of A q = lambda B q.
 *
 * Uses six global rigid-motion fields and 18 quadratic displacement fields.
 * All fields are sampled at independent master DOFs, so the trial space respects
 * the homogeneous constraints. Nodal rotations use half the displacement curl.
 *
 * A must be positive semidefinite; B may be indefinite (e.g. for buckling).
 * Eigen solves only dense projected problems of dimension at most 24.
 * The result is a Ritz estimate, not a guaranteed accurate shift or a complete
 * nullspace analysis. NaN signals that no usable positive estimate was found.
 */
Precision estimate_eigenvalue_ritz(
    const SparseMatrix&                      A,
    const SparseMatrix&                      B,
    const model::Field&                      positions,
    const constraint::ConstraintTransformer& transformer) {

    const Precision invalid   = std::numeric_limits<Precision>::quiet_NaN();
    const Precision tolerance = Precision(100) * std::numeric_limits<Precision>::epsilon();
    const auto&     map       = transformer.reduced_map();

    // ---------------------------------------------------------------------
    // (1) Validate the reduced system and access prepared master locations
    // ---------------------------------------------------------------------
    // The matrices are already constraint-reduced. The transformer supplies
    // the physical interpretation of each independent coordinate.
    logging::error(A.rows() == A.cols(),
        "Ritz estimator requires a square stiffness matrix");
    logging::error(B.rows() == A.rows() && B.cols() == A.cols(),
        "Ritz estimator requires matching operator dimensions");
    logging::error(map.n_master() == static_cast<Index>(A.rows()),
        "Ritz estimator: master count does not match reduced system");
    logging::error(map.full_size == static_cast<Index>(map.T.rows()) && map.T.cols() == A.rows(),
        "Ritz estimator: invalid constraint transformation dimensions");
    logging::error(transformer.has_dof_locations(),
        "Ritz estimator requires initialized transformer DOF locations");
    logging::error(positions.components >= 3,
        "Ritz estimator requires three-dimensional positions");

    if (A.rows() == 0) return invalid;

    // The bounding box includes all active nodes, including nodes whose DOFs
    // are fully constrained. Unused node rows do not influence the span.
    Vec3 minimum = Vec3::Constant(std::numeric_limits<Precision>::infinity());
    Vec3 maximum = -minimum;

    for (Index node : transformer.active_nodes()) {
        logging::error(node < positions.rows,
            "Ritz estimator: active node outside position field: ", node);

        const Vec3 x(positions(node, 0), positions(node, 1), positions(node, 2));
        if (!x.allFinite()) return invalid;
        minimum = minimum.cwiseMin(x);
        maximum = maximum.cwiseMax(x);
    }

    // ---------------------------------------------------------------------
    // (2) Sample six linear and 18 quadratic displacement fields
    // ---------------------------------------------------------------------
    // First six columns:
    //     translations u = e_i, rotations u = e_i x (X - center).
    // Remaining columns: for f = x²,y²,z²,xy,xz,yz and each direction e_i,
    //     u = L f e_i, theta = (L/2) grad_X(f) x e_i.
    // Coordinates are dimensionless inside f, but gradients are physical.
    // Every sampled master vector defines a valid reduced displacement q;
    // dependent slave values are recovered implicitly by the constraints.
    const Vec3      center    = minimum / 2 + maximum / 2;
    const Vec3      half_span = maximum / 2 - minimum / 2;
    const Precision length    = 2 * half_span.maxCoeff();
    if (!std::isfinite(length) || length <= 0) return invalid;

    Vec3 inverse_span = Vec3::Zero();
    for (int axis = 0; axis < 3; ++axis) {
        if (half_span[axis] > 0) inverse_span[axis] = Precision(1) / half_span[axis];
    }

    DynamicMatrix V = DynamicMatrix::Zero(A.rows(), 24);
    for (Index master = 0; master < map.n_master(); ++master) {
        const auto [node, dof] = transformer.master_location(master);
        logging::error(dof >= 0 && dof < 6,
            "Ritz estimator requires translational/rotational DOFs");
        const Vec3  r(positions(node, 0) - center[0],
                      positions(node, 1) - center[1],
                      positions(node, 2) - center[2]);
        const Vec3 p = r.cwiseProduct(inverse_span);

        // A global rigid-body rotation has both translational and rotational
        // DOFs. Unlike the old implementation, we need not pre-classify these
        // columns as rigid: if supports make them flexible, Ritz keeps them.
        if (dof < 3) {
            V(master, dof) = 1;
            for (int axis = 0; axis < 3; ++axis) {
                V(master, 3 + axis) = Vec3::Unit(axis).cross(r)[dof];
            }
        } else {
            V(master, dof) = 1;
        }

        const Precision x = p[0], y = p[1], z = p[2];
        const Precision ix = inverse_span[0], iy = inverse_span[1], iz = inverse_span[2];
        const std::array<Precision, 6> values{
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

        for (int polynomial = 0; polynomial < 6; ++polynomial) {
            for (int direction = 0; direction < 3; ++direction) {
                const int column = 6 + 3 * polynomial + direction;
                if (dof < 3) {
                    if (dof == direction) V(master, column) = length * values[polynomial];
                } else {
                    V(master, column) = Precision(0.5) * length
                        * gradients[polynomial].cross(Vec3::Unit(direction))[dof - 3];
                }
            }
        }
    }
    if (!V.allFinite()) return invalid;

    // ---------------------------------------------------------------------
    // (3) Build an independent Euclidean-orthonormal trial basis Q
    // ---------------------------------------------------------------------
    // Normalize each candidate so that the Gram matrix V^T V can identify
    // dependent fields without being distorted by their arbitrary amplitudes.
    // Its eigendecomposition replaces modified Gram-Schmidt on full-size DOFs.
    for (Eigen::Index j = 0; j < V.cols(); ++j) {
        const Precision norm = V.col(j).stableNorm();
        if (!std::isfinite(norm)) return invalid;
        if (norm > 0) V.col(j) /= norm;
    }

    Eigen::SelfAdjointEigenSolver<DynamicMatrix> gram_eig(symmetric(V.transpose() * V));
    if (gram_eig.info() != Eigen::Success) return invalid;

    const Precision gram_scale = gram_eig.eigenvalues().cwiseAbs().maxCoeff();
    if (!(gram_scale > 0)) return invalid;

    Eigen::Index count = 0;
    for (Eigen::Index i = 0; i < gram_eig.eigenvalues().size(); ++i) {
        if (gram_eig.eigenvalues()[i] > tolerance * gram_scale) ++count;
        if (gram_eig.eigenvalues()[i] < -tolerance * gram_scale) return invalid;
    }
    if (count == 0) return invalid;

    DynamicMatrix transform(V.cols(), count);
    Eigen::Index column = 0;
    for (Eigen::Index i = 0; i < gram_eig.eigenvalues().size(); ++i) {
        const Precision value = gram_eig.eigenvalues()[i];
        if (value <= tolerance * gram_scale) continue;
        transform.col(column++) = gram_eig.eigenvectors().col(i) / std::sqrt(value);
    }
    const DynamicMatrix Q = V * transform;

    // ---------------------------------------------------------------------
    // (4) Reduce the full operators and separate stiff and zero-energy modes
    // ---------------------------------------------------------------------
    // A = Q^T A Q must be positive semidefinite. Eigen's small symmetric
    // eigensolver splits its nullspace Z from positive-energy directions F.
    // Scale F such that F^T A F = I, without forming an inverse of the FE K.
    const DynamicMatrix Ar = symmetric(Q.transpose() * (A * Q));
    const DynamicMatrix Br = symmetric(Q.transpose() * (B * Q));
    if (!Ar.allFinite() || !Br.allFinite()) return invalid;

    Eigen::SelfAdjointEigenSolver<DynamicMatrix> stiffness_eig(Ar);
    if (stiffness_eig.info() != Eigen::Success) return invalid;

    const Precision stiffness_scale = stiffness_eig.eigenvalues().cwiseAbs().maxCoeff();
    if (!(stiffness_scale > 0)) return invalid;

    DynamicMatrix F(Q.cols(), Q.cols());
    DynamicMatrix Z(Q.cols(), Q.cols());
    Eigen::Index flexible_count = 0;
    Eigen::Index zero_count     = 0;

    for (Eigen::Index i = 0; i < Q.cols(); ++i) {
        const Precision value = stiffness_eig.eigenvalues()[i];
        if (value < -tolerance * stiffness_scale) return invalid;
        if (value > tolerance * stiffness_scale) {
            F.col(flexible_count++) = stiffness_eig.eigenvectors().col(i) / std::sqrt(value);
        } else {
            Z.col(zero_count++) = stiffness_eig.eigenvectors().col(i);
        }
    }
    if (flexible_count == 0) return invalid;
    F.conservativeResize(Eigen::NoChange, flexible_count);
    Z.conservativeResize(Eigen::NoChange, zero_count);

    // ---------------------------------------------------------------------
    // (5) Eliminate the zero-stiffness directions in the small problem
    // ---------------------------------------------------------------------
    // With q = F b + Z a, the reciprocal eigenproblem B q = mu A q reads
    //     G a + H b = 0,
    //     H^T a + C b = mu b,
    // where G = Z^T B Z, H = Z^T B F, C = F^T B F.
    // If G is invertible, eliminate a to obtain C - H^T G^-1 H.
    // If G is singular (possible for geometric stiffness), its nullspace
    // instead imposes H_null b = 0. A small SVD enforces this constraint.
    DynamicMatrix effective_B = symmetric(F.transpose() * Br * F);

    if (zero_count > 0) {
        const DynamicMatrix G = symmetric(Z.transpose() * Br * Z);
        const DynamicMatrix H = Z.transpose() * Br * F;

        Eigen::SelfAdjointEigenSolver<DynamicMatrix> rigid_eig(G);
        if (rigid_eig.info() != Eigen::Success) return invalid;

        // Compare against the entire B operator: a true zero of G can
        // acquire a roundoff-sized eigenvalue even if G has no nonzero scale.
        const Precision metric_scale = Br.norm();
        DynamicMatrix constraints(zero_count, flexible_count);
        Eigen::Index constraint_count = 0;

        for (Eigen::Index i = 0; i < zero_count; ++i) {
            const Precision value = rigid_eig.eigenvalues()[i];
            const Eigen::Matrix<Precision, 1, Eigen::Dynamic> h = rigid_eig.eigenvectors().col(i).transpose() * H;

            if (std::abs(value) > tolerance * metric_scale) {
                effective_B.noalias() -= h.transpose() * (h / value);
            } else {
                constraints.row(constraint_count++) = h;
            }
        }

        // Null-metric rigid directions can couple to flexible fields even
        // though their own B energy vanishes. The admissible flexible space
        // is the kernel of these remaining coupling constraints.
        if (constraint_count > 0) {
            constraints.conservativeResize(constraint_count, Eigen::NoChange);
            Eigen::JacobiSVD<DynamicMatrix> svd(constraints, Eigen::ComputeFullV);
            if (svd.info() != Eigen::Success) return invalid;

            // Use the underlying B scale as well, so roundoff in an
            // identically zero coupling does not become a constraint.
            const Precision singular_scale = Br.norm() * F.norm();
            Eigen::Index rank = 0;
            for (Eigen::Index i = 0; i < svd.singularValues().size(); ++i) {
                if (svd.singularValues()[i] > tolerance * singular_scale) ++rank;
            }
            if (rank == flexible_count) return invalid;

            const DynamicMatrix N = svd.matrixV().rightCols(flexible_count - rank);
            effective_B = symmetric(N.transpose() * effective_B * N);
        }
    }

    // ---------------------------------------------------------------------
    // (6) The largest positive reciprocal eigenvalue gives the estimate
    // ---------------------------------------------------------------------
    // In the remaining positive-stiffness space A = I, so only the ordinary
    // symmetric eigenproblem effective_B b = mu b remains. The smallest
    // positive Ritz eigenvalue of the original problem is lambda = 1 / mu_max.
    if (!effective_B.allFinite()) return invalid;
    Eigen::SelfAdjointEigenSolver<DynamicMatrix> eig(symmetric(effective_B), Eigen::EigenvaluesOnly);
    if (eig.info() != Eigen::Success) return invalid;

    const Precision reciprocal = eig.eigenvalues().maxCoeff();
    if (!std::isfinite(reciprocal) || reciprocal <= 0) return invalid;

    const Precision estimate = Precision(1) / reciprocal;
    return std::isfinite(estimate) && estimate > 0 ? estimate : invalid;
}

} // namespace mattools
} // namespace fem
