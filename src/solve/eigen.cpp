/******************************************************************************
 * @file eigen.cpp
 * @brief Public dispatchers for sparse symmetric eigenvalue solvers.
 *
 * @details This file exposes two overloads of `eigs(...)`:
 *          - Simple problem:       A x = λ x
 *          - Generalized problem:  A x = λ B x
 *
 *          It validates inputs, ensures compressed storage, and dispatches to
 *          the appropriate implementation unit based on the selected mode:
 *          - Simple:      Regular, ShiftInvert
 *          - Generalized: ShiftInvert, Buckling, Cayley
 *
 *          Any unsupported mode is rejected with a clear error.
 *
 * @author Created by Finn Eggers (c) <finn.eggers@rwth-aachen.de>
 *         All rights reserved.
 * @date   Created on 19.09.2025
 ******************************************************************************/

#include "eigen_internal.h"

namespace fem::solver {

// Simple: A x = λ x
std::vector<EigenValueVectorPair>
eigs(SolverDevice device, const SparseMatrix& A, int k, const EigenOpts& opts)
{
    logging::error(A.rows() == A.cols(), "Matrix must be square");
    logging::error(k > 0, "Requested k must be > 0");
    logging::error(k < A.rows(), "Requesting k >= N is not supported in sparse partial mode");

    if (device != CPU) {
        logging::info(true, "Eigen: GPU path not implemented; falling back to CPU");
    }

    detail::ensure_compressed(A);

    switch (opts.mode) {
        case EigenMode::Regular:
            return detail::eigs_simple_regular(A, k, opts);
        case EigenMode::ShiftInvert:
            return detail::eigs_simple_shiftinvert(A, k, opts);
        default:
            logging::error(false, "Simple eigenproblem supports only Regular or ShiftInvert");
            return {};
    }
}

// Generalized: A x = λ B x
std::vector<EigenValueVectorPair>
eigs(SolverDevice device, const SparseMatrix& A, const SparseMatrix& B, int k, const EigenOpts& opts)
{
    logging::error(A.rows() == A.cols(), "A must be square");
    logging::error(B.rows() == B.cols(), "B must be square");
    logging::error(A.rows() == B.rows(), "A and B must have the same size");
    logging::error(k > 0, "Requested k must be > 0");
    logging::error(k < A.rows(), "Requesting k >= N is not supported in sparse partial mode");

    if (device != CPU) {
        logging::info(true, "Eigen (generalized): GPU path not implemented; falling back to CPU");
    }

    detail::ensure_compressed(A);
    detail::ensure_compressed(B);

    switch (opts.mode) {
        case EigenMode::ShiftInvert:
            return detail::eigs_general_shiftinvert(A, B, k, opts);
        case EigenMode::Buckling:
            return detail::eigs_general_buckling(A, B, k, opts);
        case EigenMode::Cayley:
            return detail::eigs_general_cayley(A, B, k, opts);
        default:
            logging::error(false, "Generalized eigenproblem supports only ShiftInvert, Buckling, or Cayley");
            return {};
    }
}

} // namespace fem::solver
