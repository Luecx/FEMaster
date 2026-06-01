#pragma once

#include "solve_newmark.h"

namespace fem::solver::detail {

NewmarkResult
newmark_linear_gpu(const SparseMatrix& M,
                   const SparseMatrix& C,
                   const SparseMatrix& K,
                   const NewmarkIC& ic,
                   const NewmarkOpts& opts,
                   ForceFn f_of_t);

} // namespace fem::solver::detail
