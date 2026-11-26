/**
 * @file extract_connected_rows.cpp
 * @brief Implements extraction of rows coupled to seed rows through shared columns.
 *
 * A breadth-first search alternates between rows and columns of the matrix.
 * Each row and each column is visited at most once, so the search is linear
 * in the number of stored entries.
 *
 * @see extract_connected_rows
 * @author Finn Eggers
 * @date 25.11.2026
 */

#include "extract_connected_rows.h"

#include "../core/logging.h"

#include <algorithm>

namespace fem {
namespace mattools {

/**
 * Extracts the row block connected to the given seed rows.
 *
 * A row-major copy provides the columns of each row, the column-major source
 * provides the rows of each column. The selected rows act as the search queue.
 * Explicitly stored zeros do not couple rows but are copied for selected rows.
 *
 * @param matrix Source matrix, unchanged.
 * @param seeds Start rows of the search.
 * @param rows Receives the original indices of the selected rows, ascending.
 * @return Selected rows with all original columns.
 */
SparseMatrix extract_connected_rows(const SparseMatrix&       matrix,
                                    const std::vector<Index>& seeds,
                                    std::vector<Index>&       rows) {
    using RowMajorMatrix = Eigen::SparseMatrix<Precision, Eigen::RowMajor>;

    const Index n_rows = matrix.rows();
    const Index n_cols = matrix.cols();

    std::vector<bool> row_visited(n_rows, false);
    std::vector<bool> col_visited(n_cols, false);

    rows.clear();
    for (const Index seed : seeds) {
        logging::error(seed >= 0 && seed < n_rows,
            "extract_connected_rows: seed row ", seed, " out of range");
        if (row_visited[seed]) continue;
        row_visited[seed] = true;
        rows.push_back(seed);
    }

    const RowMajorMatrix by_row = matrix;
    for (Index next = 0; next < static_cast<Index>(rows.size()); ++next) {
        for (RowMajorMatrix::InnerIterator row_entry(by_row, rows[next]); row_entry; ++row_entry) {
            const Index col = row_entry.col();
            if (row_entry.value() == Precision(0) || col_visited[col]) continue;
            col_visited[col] = true;

            for (SparseMatrix::InnerIterator col_entry(matrix, col); col_entry; ++col_entry) {
                const Index row = col_entry.row();
                if (col_entry.value() == Precision(0) || row_visited[row]) continue;
                row_visited[row] = true;
                rows.push_back(row);
            }
        }
    }
    std::sort(rows.begin(), rows.end());

    // Map selected rows to their local position and copy their entries.
    std::vector<Index> local_row(n_rows, 0);
    for (Index i = 0; i < static_cast<Index>(rows.size()); ++i) {
        local_row[rows[i]] = i;
    }

    TripletList triplets;
    for (Eigen::Index col = 0; col < matrix.outerSize(); ++col) {
        for (SparseMatrix::InnerIterator entry(matrix, col); entry; ++entry) {
            const Index row = entry.row();
            if (!row_visited[row]) continue;
            triplets.emplace_back(static_cast<SparseMatrix::StorageIndex>(local_row[row]),
                                  static_cast<SparseMatrix::StorageIndex>(col),
                                  entry.value());
        }
    }

    SparseMatrix result(static_cast<Eigen::Index>(rows.size()), matrix.cols());
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

} // namespace mattools
} // namespace fem
