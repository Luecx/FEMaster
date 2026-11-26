#pragma once

#include "solve_eigval_common.h"

namespace fem::solver::detail {

std::vector<EigvalPair>
eigval_simple_cpu(const SparseMatrix& A, Index k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_general_cpu(const SparseMatrix& A, const SparseMatrix& B, Index k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_general_cpu(const SparseMatrix& A, const SparseMatrix& B, Precision min_val, Precision max_val,
                   const EigvalOpts& opts, bool* interval_unresolved);

std::vector<EigvalPair>
eigval_simple_regular_cpu(const SparseMatrix& A, Index k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_simple_shiftinvert_cpu(const SparseMatrix& A, Index k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_general_shiftinvert_cpu(const SparseMatrix& A, const SparseMatrix& B, Index k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_general_buckling_cpu(const SparseMatrix& A, const SparseMatrix& B, Index k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_general_cayley_cpu(const SparseMatrix& A, const SparseMatrix& B, Index k, const EigvalOpts& opts);

} // namespace fem::solver::detail
