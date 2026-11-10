/**
 * @file projection.h
 * @brief Builds the mass-orthogonal complement of a supplied modal subspace.
 *
 * Given a matrix R whose columns span the modes to be removed and a symmetric
 * mass matrix M, the projector is
 *
 *     P = I - R (R^T M R)^(-1) R^T M.
 *
 * It satisfies P R = 0, P^2 = P and P^T M = M P, provided the normalized Gram
 * matrix R^T M R is positive definite. In particular, the columns of R do
 * not need to be mass-normalized or mutually orthogonal.
 *
 * The input directions are independently mass-normalized before the rank test.\n * The returned matrix is dense and has the same dimension as M. For large FE
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
inline DynamicMatrix build_projection_matrix(const DynamicMatrix& modes,
                                           const SparseMatrix& mass) {
    // Validate common coordinates before constructing the normalized Gram matrix.
    logging::error(mass.rows() == mass.cols(),
        "Projection matrix requires a square mass matrix");
    logging::error(modes.rows() == mass.rows(),
        "Projection matrix: mode and mass dimensions do not match");
    logging::error(modes.allFinite(),
        "Projection matrix requires finite modal vectors");

    const Eigen::Index n = mass.rows();
    DynamicMatrix projector = DynamicMatrix::Identity(n, n);

    // An empty modal basis leaves the entire displacement space unchanged.
    if (modes.cols() == 0) {
        return projector;
    }

    // Scale each direction independently before forming the Gram matrix.
    // Otherwise the rank test depends on arbitrary input mode amplitudes.
    // Initial Euclidean scaling avoids overflow in the mass quadratic form.
    DynamicMatrix normalized = modes;
    for (Eigen::Index col = 0; col < normalized.cols(); ++col) {
        const Precision scale = normalized.col(col).cwiseAbs().maxCoeff();
        logging::error(std::isfinite(scale) && scale > Precision(0),
            "Projection matrix requires nonzero modal directions");

        normalized.col(col) /= scale;
    }

    DynamicMatrix weighted = mass * normalized;
    logging::error(weighted.allFinite(),
        "Projection matrix: nonfinite mass-weighted directions");

    for (Eigen::Index col = 0; col < normalized.cols(); ++col) {
        const Precision modal_mass = normalized.col(col).dot(weighted.col(col));
        logging::error(std::isfinite(modal_mass) && modal_mass > Precision(0),
            "Projection matrix requires directions with positive modal mass");

        const Precision inverse_norm = Precision(1) / std::sqrt(modal_mass);
        normalized.col(col) *= inverse_norm;
        weighted.col(col) *= inverse_norm;
    }

    // The normalized Gram matrix has unit diagonal. Its eigenvalues therefore
    // measure dependence among directions independently of column scaling.
    const DynamicMatrix gram = normalized.transpose() * weighted;
    logging::error(gram.allFinite(),
        "Projection matrix: nonfinite normalized Gram matrix");

    Eigen::SelfAdjointEigenSolver<DynamicMatrix> eigen(gram);
    logging::error(eigen.info() == Eigen::Success,
        "Projection matrix: Gram eigendecomposition failed");

    const Precision maximum = eigen.eigenvalues().cwiseAbs().maxCoeff();
    const Precision tolerance = Precision(100) * std::numeric_limits<Precision>::epsilon()
                              * static_cast<Precision>(std::max<Eigen::Index>(1, modes.cols()))
                              * maximum;

    logging::error(maximum > Precision(0)
                   && (eigen.eigenvalues().array() > tolerance).all(),
        "Projection matrix requires linearly independent positive-mass directions");

    // Invert only the small Gram matrix; reconstruct the projector using the
    // normalized basis to avoid undoing the numerical conditioning improvement.
    const DynamicMatrix gram_inverse = eigen.eigenvectors()
        * eigen.eigenvalues().cwiseInverse().asDiagonal()
        * eigen.eigenvectors().transpose();

    // P is mass-orthogonal rather than necessarily Euclidean-symmetric.
    projector.noalias() -= normalized * gram_inverse * weighted.transpose();

    return projector;
}

} // namespace mattools
} // namespace fem
