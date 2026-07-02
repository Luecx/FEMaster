#include "regularise_stiffness.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace fem {

int regularise_stiffness(SparseMatrix& matrix, Precision alpha) {
    if (alpha <= Precision(0) || matrix.rows() == 0) {
        return 0;
    }

    Precision     diagonal_sum   = Precision(0);
    Index         diagonal_count = 0;
    Precision     entry_sum      = Precision(0);
    Index         entry_count    = 0;
    DynamicVector row_norm       = DynamicVector::Zero(matrix.rows());

    for (int col = 0; col < matrix.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator it(matrix, col); it; ++it) {
            const Precision value_abs = std::abs(it.value());

            row_norm(it.row()) += value_abs;

            if (value_abs > Precision(0)) {
                entry_sum += value_abs;
                ++entry_count;
            }

            if (it.row() == it.col() && value_abs > Precision(0)) {
                diagonal_sum += value_abs;
                ++diagonal_count;
            }
        }
    }

    Precision stiffness_scale = Precision(0);

    if (diagonal_count > 0) {
        stiffness_scale = diagonal_sum / static_cast<Precision>(diagonal_count);
    } else if (entry_count > 0) {
        stiffness_scale = entry_sum / static_cast<Precision>(entry_count);
    }

    if (stiffness_scale <= Precision(0) || !std::isfinite(stiffness_scale)) {
        return 0;
    }

    const Precision min_row_stiffness = std::max(
        alpha * stiffness_scale,
        std::numeric_limits<Precision>::epsilon()
    );

    TripletList additions;

    for (int row = 0; row < matrix.rows(); ++row) {
        if (row_norm(row) < min_row_stiffness) {
            additions.emplace_back(row, row, min_row_stiffness - row_norm(row));
        }
    }

    if (additions.empty()) {
        return 0;
    }

    SparseMatrix regularisation(matrix.rows(), matrix.cols());
    regularisation.setFromTriplets(additions.begin(), additions.end());

    matrix += regularisation;
    matrix.makeCompressed();

    return static_cast<int>(additions.size());
}

} // namespace fem
