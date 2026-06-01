#include "solve_newmark_gpu.h"

#include "../../core/logging.h"

namespace fem::solver::detail {

NewmarkResult
newmark_linear_gpu(const SparseMatrix&,
                   const SparseMatrix&,
                   const SparseMatrix&,
                   const NewmarkIC&,
                   const NewmarkOpts&,
                   ForceFn)
{
    logging::error(false, "Newmark: GPU path is not implemented yet.");
    return {};
}

} // namespace fem::solver::detail
