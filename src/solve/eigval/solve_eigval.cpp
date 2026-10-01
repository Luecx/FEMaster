#include "solve_eigval_cpu.h"
#include "solve_eigval_gpu.h"

#include <cmath>

namespace fem::solver {

std::vector<EigvalPair>
eigvals(SolverDevice device, const SparseMatrix& A, int k, const EigvalOpts& opts) {
    logging::error(A.rows() == A.cols(), "Matrix must be square");
    logging::error(k > 0, "Requested k must be > 0");
    logging::error(k < A.rows(), "Requesting k >= N is not supported in sparse partial mode");

    detail::ensure_compressed(A);

    return (device == GPU)
        ? detail::eigval_simple_gpu(A, k, opts)
        : detail::eigval_simple_cpu(A, k, opts);
}

std::vector<EigvalPair>
eigvals(SolverDevice device, const SparseMatrix& A, const SparseMatrix& B, int k, const EigvalOpts& opts) {
    logging::error(A.rows() == A.cols(), "A must be square");
    logging::error(B.rows() == B.cols(), "B must be square");
    logging::error(A.rows() == B.rows(), "A and B must have the same size");
    logging::error(k > 0, "Requested k must be > 0");
    logging::error(k < A.rows(), "Requesting k >= N is not supported in sparse partial mode");

    detail::ensure_compressed(A);
    detail::ensure_compressed(B);

    return (device == GPU)
        ? detail::eigval_general_gpu(A, B, k, opts)
        : detail::eigval_general_cpu(A, B, k, opts);
}

std::vector<EigvalPair>
eigvals(SolverDevice device, const SparseMatrix& A, const SparseMatrix& B, Precision min_val, Precision max_val, const EigvalOpts& opts) {
    logging::error(A.rows() == A.cols(), "A must be square");
    logging::error(B.rows() == B.cols(), "B must be square");
    logging::error(A.rows() == B.rows(), "A and B must have the same size");
    logging::error(std::isfinite(min_val) && std::isfinite(max_val) && min_val < max_val,
        "Requested eigenvalue bounds must be finite and satisfy min_val < max_val");

    detail::ensure_compressed(A);
    detail::ensure_compressed(B);

    return (device == GPU)
        ? detail::eigval_general_gpu(A, B, min_val, max_val, opts)
        : detail::eigval_general_cpu(A, B, min_val, max_val, opts);
}


} // namespace fem::solver
