/**
 * @file solve_eigval_cpu.cpp
 * @brief Dispatches CPU eigenvalue solves and extracts generalized modal intervals.
 *
 * Count-based solves are forwarded to the regular, shift-invert, buckling or
 * Cayley implementation. Interval extraction uses generalized shift-invert and
 * retains eigenpairs strictly inside the supplied eigenvalue bounds. With MKL,
 * the inertia of A - max_val B determines how many eigenvalues lie below the
 * upper bound; other CPU backends enlarge the requested count iteratively.
 *
 * Matrix assembly, constraint reduction and physical frequency scaling belong
 * to the caller. Interval extraction preserves the caller's spectral shift so
 * free-body stiffness matrices need not be factorized at zero.
 *
 * @see EigvalOpts
 * @author Finn Eggers
 * @date 27.10.2026
 */

#include "solve_eigval_cpu.h"
#include <algorithm>

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

std::vector<EigvalPair>
eigval_general_cpu(const SparseMatrix& A, const SparseMatrix& B, Precision min_val, Precision max_val, const EigvalOpts& opts) {
    logging::error(opts.mode == EigvalMode::ShiftInvert,
        "Generalized eigval interval supports only ShiftInvert for eigenfrequency analysis");
    logging::error(A.rows() > 1,
        "Generalized eigval interval requires N > 1 in sparse partial mode");
    logging::error(min_val < max_val,
        "Generalized eigval interval requires min_val < max_val");

    // Keep the caller's shift: modal analysis places it below the nonnegative
    // spectrum. Forcing zero here would factor a singular free-body stiffness
    // and undo the modal frequency scaling and shift selection.
    EigvalOpts opts_new = opts;
    opts_new.sort  = EigvalOpts::Sort::LargestMagn;

    const int max_attempt = static_cast<int>(A.rows()) - 1;
#ifdef USE_MKL
    // compute eigenvalues below min and below max
    Eigen::PardisoLDLT<SparseMatrix> solver{};
    SparseMatrix D = A - max_val * B;
    D.makeCompressed();
    solver.compute(D);
    logging::error(solver.info() == Eigen::Success,
        "Generalized eigval interval: shifted stiffness factorization failed");
    const int num_below_max = solver.pardisoParameterArray()[22];

    if (num_below_max == 0) {
        return {};
    }
    logging::error(num_below_max > 0 && num_below_max <= max_attempt,
        "Generalized eigval interval requires k < N in sparse partial mode");

    const int k = num_below_max + std::min(num_below_max / 5, max_attempt - num_below_max);
    std::vector<EigvalPair> eig_values_raw = eigval_general_cpu(A, B, k, opts_new);
#else
    std::vector<EigvalPair> eig_values_raw;
    const int min_attempt = std::min(10, max_attempt);
    for (int i = min_attempt; i <= max_attempt; i += std::min(i, max_attempt - i)) {
        eig_values_raw = eigval_general_cpu(A, B, i, opts_new);
        // check if the maximum eigenvalue is above max, if so, we dont need to look further
        const bool interval_covered = std::any_of(eig_values_raw.begin(), eig_values_raw.end(),
            [&](const EigvalPair& e) { return e.value >= max_val; });
        if (interval_covered) {
            break;
        }

        // break if upper limit is reached
        logging::error(i < max_attempt,
            "Generalized eigval interval cannot establish complete coverage with k < N");
    }
#endif

    // filter all that are not in range
    std::vector<EigvalPair> eig_values;
    for (const auto& e : eig_values_raw) {
        if (e.value < max_val && e.value > min_val) eig_values.push_back(e);
    }
    return eig_values;

}


} // namespace fem::solver::detail
