#pragma once

#include "solve_eigval_common.h"

namespace fem::solver::detail {

std::vector<EigvalPair>
eigval_simple_gpu(const SparseMatrix& A, int k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_general_gpu(const SparseMatrix& A, const SparseMatrix& B, int k, const EigvalOpts& opts);

} // namespace fem::solver::detail
