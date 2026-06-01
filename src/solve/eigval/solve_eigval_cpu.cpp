#include "solve_eigval_cpu.h"

namespace fem::solver::detail {

std::vector<EigvalPair>
eigval_simple_cpu(const SparseMatrix& A, int k, const EigvalOpts& opts) {
    switch (opts.mode) {
        case EigvalMode::Regular:
            return eigval_simple_regular_cpu(A, k, opts);
        case EigvalMode::ShiftInvert:
            return eigval_simple_shiftinvert_cpu(A, k, opts);
        default:
            logging::error(false, "Simple eigval problem supports only Regular or ShiftInvert");
            return {};
    }
}

std::vector<EigvalPair>
eigval_general_cpu(const SparseMatrix& A, const SparseMatrix& B, int k, const EigvalOpts& opts) {
    switch (opts.mode) {
        case EigvalMode::ShiftInvert:
            return eigval_general_shiftinvert_cpu(A, B, k, opts);
        case EigvalMode::Buckling:
            return eigval_general_buckling_cpu(A, B, k, opts);
        case EigvalMode::Cayley:
            return eigval_general_cayley_cpu(A, B, k, opts);
        default:
            logging::error(false, "Generalized eigval problem supports only ShiftInvert, Buckling, or Cayley");
            return {};
    }
}

} // namespace fem::solver::detail
