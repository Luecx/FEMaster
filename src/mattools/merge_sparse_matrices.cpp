/**
 * @file merge_sparse_matrices.cpp
 * @brief Implements stable column-wise reduction of sparse assembly matrices.
 *
 * A local loser-tree traversal is shared by symbolic counting and numerical
 * accumulation. The symbolic pass determines the union of stored coordinates;
 * a prefix sum allocates exact CSC storage, and the numeric pass writes sorted
 * rows and sums overlapping values. Element evaluation and partial-matrix
 * construction remain responsibilities of the assembly routine.
 *
 * @see merge_sparse_matrices
 */

#include "merge_sparse_matrices.h"

#include "../core/logging.h"

#include <limits>

namespace fem::mattools {

/**
 * Sums compressed worker-local matrices into a sorted global CSC matrix.
 *
 * Each output column is traversed twice using the same local loser-tree
 * operation. Only nonempty input columns enter its tournament. Internal nodes
 * retain the losing input, so consuming a winner replays only its path to the
 * root. Equal row keys are ordered by input-matrix index, preserving the
 * floating-point accumulation order for a fixed collection of partial matrices.
 *
 * The first traversal counts unique rows without reading numerical values.
 * A prefix sum converts the counts into output offsets. The second traversal
 * writes each row once and accumulates all its contributions, retaining explicit
 * zeros and cancellations. Empty columns have equal successive output offsets.
 * Fixed arrays bound the tournament to 128 inputs and avoid per-column heap
 * allocation. OpenMP assigns disjoint output columns to independent workers.
 *
 * @param matrices One to 128 identically sized square matrices. Inputs must be
 *                 compressed with sorted row indices and remain unchanged.
 * @return Sum of the input matrices with a compressed, sorted CSC pattern.
 */
SparseMatrix merge_sparse_matrices(const std::vector<SparseMatrix>& matrices) {
    using StorageIndex = SparseMatrix::StorageIndex;

    // Validate the input count against the fixed tournament capacity
    constexpr int max_inputs = 128;
    logging::error(!matrices.empty() && matrices.size() <= max_inputs,
        "merge_sparse_matrices: requires between 1 and 128 input matrices");

    // Prepare one structural count per output column and its terminal offset
    const int num_inputs  = static_cast<int>(matrices.size());
    const int global_size = static_cast<int>(matrices.front().rows());
    std::vector<Index> column_offsets(global_size + 1, 0);

    // Share one column traversal between counting and writing. Null output
    // pointers select structural counting without reading coefficient arrays.
    auto merge_column = [&](int col, StorageIndex* output_rows, Precision* output_values) -> Index {
        constexpr StorageIndex exhausted = std::numeric_limits<StorageIndex>::max();

        StorageIndex        positions   [max_inputs];
        StorageIndex        ends        [max_inputs];
        StorageIndex        keys        [max_inputs];
        const StorageIndex* input_rows  [max_inputs];
        const Precision*    input_values[max_inputs];
        int                 losers      [max_inputs];
        int                 winners     [max_inputs * 2];
        int                 active = 0;

        // Compact nonempty columns in input order to preserve stable ties
        for (int input = 0; input < num_inputs; ++input) {
            const auto&        matrix = matrices[input];
            const StorageIndex begin  = matrix.outerIndexPtr()[col];
            const StorageIndex end    = matrix.outerIndexPtr()[col + 1];
            if (begin < end) {
                positions[active]  = begin;
                ends[active]       = end;
                keys[active]       = matrix.innerIndexPtr()[begin];
                input_rows[active] = matrix.innerIndexPtr();
                if (output_values) {
                    input_values[active] = matrix.valuePtr();
                }
                ++active;
            }
        }
        if (active == 0) {
            return 0;
        }

        // Build the complete tournament, padding it with exhausted leaves
        int leaves = 1;
        while (leaves < active) {
            leaves *= 2;
        }
        for (int slot = 0; slot < leaves; ++slot) {
            if (slot >= active) {
                keys[slot] = exhausted;
            }
            winners[leaves + slot] = slot;
        }
        for (int node = leaves - 1; node > 0; --node) {
            const int  left      = winners[2 * node];
            const int  right     = winners[2 * node + 1];
            const bool left_wins = keys[left] < keys[right] || (keys[left] == keys[right] && left < right);
            winners[node] = left_wins ? left : right;
            losers[node]  = left_wins ? right : left;
        }
        int champion = winners[1];

        // Count or write complete groups of equal row keys
        Index count = 0;
        while (keys[champion] != exhausted) {
            const StorageIndex row = keys[champion];
            Precision          sum = Precision(0);
            do {
                if (output_values) {
                    sum += input_values[champion][positions[champion]];
                }

                // Replace the consumed head with its successor or the sentinel
                const int previous = champion;
                ++positions[previous];
                keys[previous] = positions[previous] < ends[previous]
                    ? input_rows[previous][positions[previous]] : exhausted;

                // Replay only the affected path, retaining each losing input
                int candidate = previous;
                for (int node = (leaves + previous) / 2; node > 0; node /= 2) {
                    const int opponent = losers[node];
                    if (keys[opponent] < keys[candidate] ||
                        (keys[opponent] == keys[candidate] && opponent < candidate)) {
                        losers[node] = candidate;
                        candidate    = opponent;
                    }
                }
                champion = candidate;
            } while (keys[champion] == row);

            if (output_values) {
                output_rows[count]   = row;
                output_values[count] = sum;
            }
            ++count;
        }
        return count;
    };

    // Count the structural union independently for every output column
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_inputs) schedule(static, 256)
#endif
    for (int col = 0; col < global_size; ++col) {
        column_offsets[col + 1] = merge_column(col, nullptr, nullptr);
    }

    // Convert column counts into disjoint compressed-storage ranges
    for (int col = 0; col < global_size; ++col) {
        column_offsets[col + 1] += column_offsets[col];
    }
    const Index total_nonzeros = column_offsets[global_size];

    // Allocate the exact global CSC pattern and expose its compressed arrays
    SparseMatrix result(global_size, global_size);
    result.makeCompressed();
    result.resizeNonZeros(static_cast<Eigen::Index>(total_nonzeros));

    auto* outer_indices = result.outerIndexPtr();
    auto* inner_indices = result.innerIndexPtr();
    auto* values        = result.valuePtr();

    // Accumulate each column directly into its exclusive output range.
    // Empty columns need no traversal and may have null coefficient arrays.
#ifdef _OPENMP
#pragma omp parallel for num_threads(num_inputs) schedule(static, 256)
#endif
    for (int col = 0; col < global_size; ++col) {
        const Index start = column_offsets[col];
        outer_indices[col] = static_cast<StorageIndex>(start);
        if (column_offsets[col + 1] > start) {
            merge_column(col, inner_indices + start, values + start);
        }
    }
    outer_indices[global_size] = static_cast<StorageIndex>(total_nonzeros);

    return result;
}

} // namespace fem::mattools
