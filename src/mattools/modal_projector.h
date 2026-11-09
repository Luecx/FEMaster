/**
 * @file modal_projector.h
 * @brief Builds the mass-orthogonal complement of a supplied modal subspace.
 *
 * Given a matrix R whose columns span the modes to be removed and a symmetric
 * mass matrix M, the projector is
 *
 *     P = I - R (R^T M R)^(-1) R^T M.
 *
 * It satisfies P R = 0, P^2 = P and P^T M = M P, provided the modal Gram
 * matrix R^T M R is positive definite. In particular, the columns of R do
 * not need to be mass-normalized or mutually orthogonal.
 *
 * The returned matrix is dense and has the same dimension as M. For large FE
 * models, applying P implicitly is considerably cheaper than forming P and
 * computing P^T K P explicitly. Projection does not reduce the matrix size:
 * P^T K P still has a nullspace associated with the removed directions.
 *
 * @author Finn Eggers
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

#include <Eigen/Eigenvalues>

#include "../core/logging.h"
#include "../core/types_eig.h"

namespace fem {
namespace mattools {

/**
 * @brief Construct the mass-orthogonal projector that removes given modes.
 *
 * Both M and the columns of modes must be expressed in the same coordinates,
 * e.g. after applying the constraint transformation to the FE system. Only
 * actual admissible rigid-body modes should be supplied when removing rigid
 * motion; geometric rigid directions generally violate constrained DOFs.
 *
 * This routine assumes a symmetric positive-semidefinite mass matrix and
 * requires positive modal mass for every independent column. A massless mode
 * cannot be removed by this mass-orthogonal projector.
 *
 * @param modes Full-column-rank matrix of modes to remove, one per column.
 * @param mass Symmetric mass matrix in the same coordinate system.
 * @return Dense projector P, such that P * modes is zero to roundoff.
 */
inline DynamicMatrix build_modal_projector(const DynamicMatrix& modes,
                                           const SparseMatrix& mass) {
    // Validate common coordinates before constructing the modal Gram matrix.
    logging::error(mass.rows() == mass.cols(),
        "Modal projector requires a square mass matrix");
    logging::error(modes.rows() == mass.rows(),
        "Modal projector: mode and mass dimensions do not match");
    logging::error(modes.allFinite(),
        "Modal projector requires finite modal vectors");

    const Eigen::Index n = mass.rows();
    DynamicMatrix projector = DynamicMatrix::Identity(n, n);

    // An empty modal basis leaves the entire displacement space unchanged.
    if (modes.cols() == 0) {
        return projector;
    }

    // Only a small dense matrix (number of removed modes squared) is inverted.
    const DynamicMatrix weighted_modes = mass * modes;
    const DynamicMatrix gram = modes.transpose() * weighted_modes;

    logging::error(weighted_modes.allFinite() && gram.allFinite(),
        "Modal projector: nonfinite modal mass coefficients");

    // A symmetric eigendecomposition detects massless or dependent directions
    // and avoids silently treating a singular modal Gram matrix as invertible.
    Eigen::SelfAdjointEigenSolver<DynamicMatrix> eigen(gram);
    logging::error(eigen.info() == Eigen::Success,
        "Modal projector: modal Gram eigendecomposition failed");

    const Precision maximum = eigen.eigenvalues().cwiseAbs().maxCoeff();
    const Precision tolerance = Precision(100) * std::numeric_limits<Precision>::epsilon()
                              * static_cast<Precision>(std::max<Eigen::Index>(1, modes.cols()))
                              * maximum;

    logging::error(maximum > Precision(0)
                   && (eigen.eigenvalues().array() > tolerance).all(),
        "Modal projector requires independent modes with positive modal mass");

    // Solve the small modal Gram system via its spectral decomposition.
    // This also handles unnormalized input modes and non-diagonal modal mass.
    const DynamicMatrix gram_inverse = eigen.eigenvectors()
        * eigen.eigenvalues().cwiseInverse().asDiagonal()
        * eigen.eigenvectors().transpose();

    // P is generally not Euclidean-symmetric: it is orthogonal in the M metric.
    projector.noalias() -= modes * gram_inverse * weighted_modes.transpose();

    return projector;
}

} // namespace mattools
} // namespace fem
