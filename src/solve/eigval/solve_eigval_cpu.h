#pragma once

#include "solve_eigval_common.h"

namespace fem::solver::detail {

std::vector<EigvalPair>
eigval_simple_cpu(const SparseMatrix& A, int k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_general_cpu(const SparseMatrix& A, const SparseMatrix& B, int k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_simple_regular_cpu(const SparseMatrix& A, int k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_simple_shiftinvert_cpu(const SparseMatrix& A, int k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_general_shiftinvert_cpu(const SparseMatrix& A, const SparseMatrix& B, int k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_general_buckling_cpu(const SparseMatrix& A, const SparseMatrix& B, int k, const EigvalOpts& opts);

std::vector<EigvalPair>
eigval_general_cayley_cpu(const SparseMatrix& A, const SparseMatrix& B, int k, const EigvalOpts& opts);

} // namespace fem::solver::detail
