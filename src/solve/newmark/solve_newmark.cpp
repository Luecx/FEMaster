#include "solve_newmark.h"

#include "solve_newmark_cpu.h"
#include "solve_newmark_gpu.h"

#include "../../core/logging.h"

#include <utility>

namespace fem::solver {

NewmarkResult
newmark_linear(const SparseMatrix& M,
               const SparseMatrix& C,
               const SparseMatrix& K,
               const NewmarkIC& ic,
               const NewmarkOpts& opts,
               ForceFn f_of_t,
               const NewmarkForceBasis& force_basis)
{
    logging::error(M.rows() == M.cols(), "Newmark: M must be square.");
    logging::error(C.rows() == C.cols(), "Newmark: C must be square.");
    logging::error(K.rows() == K.cols(), "Newmark: K must be square.");
    logging::error(M.rows() == C.rows() && M.rows() == K.rows(), "Newmark: matrix size mismatch.");

    if (opts.device == GPU) {
        logging::error(!force_basis.empty(), "Newmark GPU path requires force_basis.");
        return detail::newmark_linear_gpu(M, C, K, ic, opts, force_basis);
    }

    return detail::newmark_linear_cpu(M, C, K, ic, opts, std::move(f_of_t));
}

} // namespace fem::solver
