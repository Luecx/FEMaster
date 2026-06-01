#include "solve_eigval_gpu.h"
#include "solve_eigval_cpu.h"

namespace fem::solver::detail {

std::vector<EigvalPair>
eigval_simple_gpu(const SparseMatrix& A, int k, const EigvalOpts& opts) {
    logging::info(true, "Eigval: GPU path not implemented; falling back to CPU");
    return eigval_simple_cpu(A, k, opts);
}

std::vector<EigvalPair>
eigval_general_gpu(const SparseMatrix& A, const SparseMatrix& B, int k, const EigvalOpts& opts) {
    logging::info(true, "Eigval (generalized): GPU path not implemented; falling back to CPU");
    return eigval_general_cpu(A, B, k, opts);
}

} // namespace fem::solver::detail
