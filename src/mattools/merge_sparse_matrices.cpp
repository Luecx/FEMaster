/**
 * @file merge_sparse_matrices.cpp
 * @brief Reduces worker-local compressed sparse assembly matrices to one CSC matrix.
 *
 * The matrix tools layer combines already-compressed partial assembly matrices
 * without changing element evaluation, triplet insertion or input compression.
 * Two inputs are merged directly through their sorted CSC columns. Three or
 * more inputs use a dense row accumulator with a hierarchical occupancy bitmap.
 *
 * Both paths first count unique stored coordinates and then write into exactly
 * sized CSC storage. Equal coordinates are accumulated in input-matrix order,
 * explicitly stored zeros are preserved, and independent output columns are
 * processed in parallel without synchronization.
 *
 * @see merge_sparse_matrices
 * @see assemble_matrix_multithreaded
 */

#include "merge_sparse_matrices.h"

#include "../core/logging.h"
#include "../core/memory.h"
#include "../core/parallel.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#ifdef _MSC_VER
    #include <intrin.h>
#endif

namespace fem::mattools {
namespace {

using StorageIndex = SparseMatrix::StorageIndex;

/**
 * @brief Accumulates stored rows of one CSC column in a reusable bitmap tree.
 *
 * Each bit in the lowest level represents a global row index. Higher levels
 * mark nonempty 64-bit words of the level immediately below, ending in a
 * single root word. New occupancy propagates upward only when a child word
 * changes from empty to nonempty.
 *
 * Traversal from the root yields row indices in ascending order and clears
 * the visited masks. The dense values remain untouched between columns;
 * their contents are valid only for rows occupied in the current column.
 * The same structure serves symbolic counting and numerical accumulation.
 */
struct ColumnBitmap {
    // Dense row coefficients and hierarchical occupancy masks, leaves first
    std::vector<Precision>                  values;
    std::vector<std::vector<std::uint64_t>> levels;

    explicit ColumnBitmap(std::size_t rows) : values(rows) {
        std::size_t words = rows / 64 + (rows % 64 != 0);

        // Build occupancy levels until a single 64-bit root word remains
        while (words != 0) {
            levels.emplace_back(words, std::uint64_t(0));
            if (words == 1) {
                break;
            }
            words = words / 64 + (words % 64 != 0);
        }
    }

    // Return true if this row has not yet been seen in the current column
    bool mark(std::size_t row) {
        std::size_t         word = row >> 6;
        const std::uint64_t bit  = std::uint64_t(1) << (row & 63);
        std::uint64_t&      leaf = levels[0][word];

        if (leaf & bit) {
            return false;
        }

        // Mark the leaf; an already occupied word needs no propagation
        const bool was_empty = (leaf == 0);
        leaf |= bit;
        if (!was_empty) {
            return true;
        }

        // A parent changes occupancy only when its child was previously empty
        for (std::size_t level = 1; level < levels.size(); ++level) {
            const std::size_t   parent           = word >> 6;
            const std::uint64_t parent_bit       = std::uint64_t(1) << (word & 63);
            std::uint64_t&      mask             = levels[level][parent];
            const bool          parent_was_empty = (mask == 0);

            mask |= parent_bit;
            if (!parent_was_empty) {
                break;
            }
            word = parent;
        }
        return true;
    }

    // Initialize the first contribution, then sum subsequent contributions
    // in input order without resetting the dense values between columns
    void add(std::size_t row, Precision value) {
        if (mark(row)) {
            values[row] = value;
        } else {
            values[row] += value;
        }
    }

    // Visit occupied rows in ascending order and clear all visited masks.
    // The callback receives only a row index so the symbolic pass never reads
    // dense coefficient values.
    template<typename Emit>
    void flush(Emit&& emit) {
        // Walk only occupied subtrees, visiting rows in ascending order
        auto visit = [&](auto&& self, std::size_t level, std::size_t word) -> void {
            std::uint64_t remaining = levels[level][word];
            levels[level][word]      = 0;

            while (remaining != 0) {
                unsigned bit = 0;
#ifdef _MSC_VER
                unsigned long position = 0;
                _BitScanForward64(&position, remaining);
                bit = static_cast<unsigned>(position);
#else
                bit = static_cast<unsigned>(__builtin_ctzll(remaining));
#endif
                const std::size_t child = (word << 6) + bit;

                // Leaf bits represent rows; higher-level bits represent child words
                if (level == 0) {
                    emit(child);
                } else {
                    self(self, level - 1, child);
                }
                remaining &= remaining - 1;
            }
        };

        // Start at the root; an empty root has no occupied bits to visit
        if (!levels.empty()) {
            visit(visit, levels.size() - 1, 0);
        }
    }
};

/**
 * Constructs the final compressed CSC matrix from per-column unique counts.
 *
 * Converts the counts into column offsets, verifies the Eigen storage-index
 * range and exact allocation size, then allocates the inner-index and value
 * arrays once. The returned matrix owns disjoint, correctly sized ranges for
 * each column, ready for concurrent writes without reallocations.
 *
 * @param rows Global matrix dimension.
 * @param offsets Unique column counts on entry; CSC prefix offsets on return.
 * @return Compressed matrix with its final storage allocated but not populated.
 */
SparseMatrix allocate_output(Eigen::Index rows, std::vector<std::size_t>& offsets) {
    const std::size_t columns   = offsets.size() - 1;
    const std::size_t max_index = static_cast<std::size_t>(std::numeric_limits<StorageIndex>::max());
    const std::size_t maximum   = std::numeric_limits<std::size_t>::max();

    // Convert independent column counts into global CSC offsets
    bool offsets_valid = true;
    for (std::size_t col = 0; col < columns; ++col) {
        if (offsets[col + 1] > max_index - offsets[col]) {
            offsets_valid = false;
            break;
        }
        offsets[col + 1] += offsets[col];
    }

    // Determine the exact memory requirement without overflowing size_t
    const std::size_t total_nnz     = offsets_valid  ? offsets.back() : 0;
    const std::size_t entry_size    = sizeof(StorageIndex) + sizeof(Precision);
    const bool        entries_fit   = total_nnz <= maximum / entry_size;
    const bool        pointers_fit  = columns < maximum / sizeof(StorageIndex);

    const std::size_t entry_bytes   = entries_fit  ? total_nnz * entry_size              : 0;
    const std::size_t pointer_bytes = pointers_fit ? (columns + 1) * sizeof(StorageIndex) : 0;
    const bool        output_fits   = entries_fit && pointers_fit && entry_bytes <= maximum - pointer_bytes;
    const std::size_t output_bytes  = output_fits  ? entry_bytes + pointer_bytes         : 0;

    const std::size_t available = memory::available_memory_bytes();
    const std::size_t budget    = available - available / 4;

    // Reject unsupported CSC sizes and allocations before modifying output
    logging::error(offsets_valid,
        "Sparse merge: output exceeds Eigen's sparse storage index range");
    logging::error(output_fits,
        "Sparse merge: CSC output allocation size overflow");
    logging::error(output_bytes <= budget,
        "Sparse merge: insufficient memory for output CSC storage (need ",
        output_bytes, " bytes, available ", available, " bytes)");

    // Allocate exact CSC arrays and publish the column pointers
    SparseMatrix result(rows, rows);
    result.resizeNonZeros(static_cast<Eigen::Index>(total_nnz));

    StorageIndex* outer = result.outerIndexPtr();
    for (std::size_t col = 0; col <= columns; ++col) {
        outer[col] = static_cast<StorageIndex>(offsets[col]);
    }
    return result;
}

/**
 * Merges exactly two compressed CSC matrices in linear time per column.
 *
 * The symbolic pass counts the union of each pair of sorted row-index
 * sequences. After exact output allocation, the numerical pass performs a
 * two-pointer traversal, copying unmatched entries and adding equal rows.
 * Explicit zeros and the original summation order are preserved.
 *
 * @param a First compressed CSC contribution.
 * @param b Second compressed CSC contribution.
 * @return Exactly allocated compressed sum of both matrices.
 */
SparseMatrix merge_sparse_pair(const SparseMatrix& a, const SparseMatrix& b) {
    // Prepare dimensions and worker count before accessing the CSC arrays
    const Eigen::Index cols          = a.cols();
    const std::size_t  column_count  = static_cast<std::size_t>(cols);
    const int          workers       = parallel::worker_count(cols, global_config.max_threads);

    // Capture input CSC arrays once for both merge passes
    const StorageIndex* a_outer = a.outerIndexPtr();
    const StorageIndex* b_outer = b.outerIndexPtr();
    const StorageIndex* a_rows  = a.innerIndexPtr();
    const StorageIndex* b_rows  = b.innerIndexPtr();
    const Precision*    a_vals  = a.valuePtr();
    const Precision*    b_vals  = b.valuePtr();

    std::vector<std::size_t> offsets(column_count + 1, 0);

    // Symbolic pass: count the union of the two sorted row-index sequences
    parallel::for_index(cols, workers, [&](Eigen::Index col, int /*worker*/) {
        StorageIndex       ia    = a_outer[col];
        StorageIndex       ib    = b_outer[col];
        const StorageIndex ea    = a_outer[col + 1];
        const StorageIndex eb    = b_outer[col + 1];
        std::size_t        count = 0;

        // Count matching rows once, then append any unmatched suffix
        while (ia < ea && ib < eb) {
            if (a_rows[ia] <= b_rows[ib]) {
                if (a_rows[ia] == b_rows[ib]) {
                    ++ib;
                }
                ++ia;
            } else {
                ++ib;
            }
            ++count;
        }
        count += static_cast<std::size_t>(ea - ia) + static_cast<std::size_t>(eb - ib);
        offsets[static_cast<std::size_t>(col) + 1] = count;
    });

    // Allocate the final CSC pattern and values from the exact column counts
    SparseMatrix result = allocate_output(a.rows(), offsets);
    StorageIndex* output_rows = result.innerIndexPtr();
    Precision*    output_vals = result.valuePtr();

    // Numerical pass: write sorted rows and sum duplicates directly in CSC
    parallel::for_index(cols, workers, [&](Eigen::Index col, int /*worker*/) {
        StorageIndex       ia  = a_outer[col];
        StorageIndex       ib  = b_outer[col];
        const StorageIndex ea  = a_outer[col + 1];
        const StorageIndex eb  = b_outer[col + 1];
        std::size_t        pos = offsets[static_cast<std::size_t>(col)];

        // Merge matching rows and preserve the order of input contributions
        while (ia < ea && ib < eb) {
            if (a_rows[ia] < b_rows[ib]) {
                output_rows[pos] = a_rows[ia];
                output_vals[pos] = a_vals[ia];
                ++ia;
            } else if (b_rows[ib] < a_rows[ia]) {
                output_rows[pos] = b_rows[ib];
                output_vals[pos] = b_vals[ib];
                ++ib;
            } else {
                output_rows[pos] = a_rows[ia];
                output_vals[pos] = a_vals[ia] + b_vals[ib];
                ++ia;
                ++ib;
            }
            ++pos;
        }

        // Copy the remaining sorted suffix of either input, if present
        while (ia < ea) {
            output_rows[pos] = a_rows[ia];
            output_vals[pos] = a_vals[ia];
            ++ia;
            ++pos;
        }
        while (ib < eb) {
            output_rows[pos] = b_rows[ib];
            output_vals[pos] = b_vals[ib];
            ++ib;
            ++pos;
        }
    });

    return result;
}

/**
 * Merges three or more compressed CSC matrices using hierarchical bitmaps.
 *
 * The routine sizes one dense-row accumulator per OpenMP worker and rejects
 * insufficient memory without reducing parallelism. During the symbolic pass, each
 * distinct stored row is marked once per column; the bitmap tree is cleared
 * without touching dense values. The numerical pass accumulates coefficients
 * in input order and traverses the tree to emit sorted, unique row indices.
 *
 * @param matrices Compressed contributions with common square dimensions.
 * @return Exactly allocated compressed sum of all contributions.
 */
SparseMatrix merge_sparse_bitmap(const std::vector<SparseMatrix>& matrices) {
    // Prepare global dimensions, index bounds and the requested parallelism
    const Eigen::Index rows         = matrices.front().rows();
    const Eigen::Index cols         = matrices.front().cols();
    const std::size_t  row_count    = static_cast<std::size_t>(rows);
    const std::size_t  column_count = static_cast<std::size_t>(cols);
    const std::size_t  maximum      = std::numeric_limits<std::size_t>::max();
    const int          workers      = parallel::worker_count(cols, global_config.max_threads);

    // Capture all input CSC arrays once, rather than retrieving pointers for
    // each input again in every column of the symbolic and numerical passes
    std::vector<const StorageIndex*> outer;
    std::vector<const StorageIndex*> inner;
    std::vector<const Precision*>    values;
    outer .reserve(matrices.size());
    inner .reserve(matrices.size());
    values.reserve(matrices.size());

    for (const SparseMatrix& matrix : matrices) {
        outer .push_back(matrix.outerIndexPtr());
        inner .push_back(matrix.innerIndexPtr());
        values.push_back(matrix.valuePtr());
    }

    // Count bitmap words level by level, ending with the single root word
    std::size_t bitmap_words = 0;
    std::size_t level_words  = row_count;
    do {
        level_words = level_words / 64 + (level_words % 64 != 0);
        bitmap_words += level_words;
    } while (level_words > 1);

    // Calculate dense-value and bitmap storage without overflowing size_t
    logging::error(row_count <= maximum / sizeof(Precision),
        "Sparse merge: dense value allocation size overflow");
    const std::size_t value_bytes = row_count * sizeof(Precision);

    logging::error(bitmap_words <= (maximum - value_bytes) / sizeof(std::uint64_t),
        "Sparse merge: bitmap workspace allocation size overflow");
    const std::size_t bytes_per_worker = value_bytes + bitmap_words * sizeof(std::uint64_t);

    // Keep a quarter of available memory outside the workspace budget
    const std::size_t available = memory::available_memory_bytes();
    const std::size_t budget    = available - available / 4;

    // Check the full worker team without multiplying allocation sizes, which
    // could overflow. Insufficient memory must not change the worker count.
    logging::error(bytes_per_worker <= budget / static_cast<std::size_t>(workers),
        "Sparse merge: insufficient memory for ", workers, " bitmap workspaces (",
        bytes_per_worker, " bytes per worker, available ", available, " bytes)");

    // Allocate one bitmap per worker and exact per-column symbolic counts
    std::vector<ColumnBitmap> bitmaps;
    bitmaps.reserve(static_cast<std::size_t>(workers));
    for (int worker = 0; worker < workers; ++worker) {
        bitmaps.emplace_back(row_count);
    }
    std::vector<std::size_t> offsets(column_count + 1, 0);

    // Symbolic pass: count each distinct stored row once in each output column
    parallel::for_index(cols, workers, [&](Eigen::Index col, int worker) {
        ColumnBitmap& bitmap = bitmaps[static_cast<std::size_t>(worker)];
        std::size_t   count  = 0;

        // Visit every input in matrix order, marking rows without reading values
        for (std::size_t input = 0; input < matrices.size(); ++input) {
            const StorageIndex  begin    = outer[input][col];
            const StorageIndex  end      = outer[input][col + 1];
            const StorageIndex* rows_ptr = inner[input];

            for (StorageIndex i = begin; i < end; ++i) {
                count += bitmap.mark(static_cast<std::size_t>(rows_ptr[i]));
            }
        }

        // Clear marked paths; symbolic counting does not inspect coefficient values
        bitmap.flush([](std::size_t) {});
        offsets[static_cast<std::size_t>(col) + 1] = count;
    });

    // Allocate the final CSC pattern and values from the exact column counts
    SparseMatrix result = allocate_output(rows, offsets);
    StorageIndex* output_rows = result.innerIndexPtr();
    Precision*    output_vals = result.valuePtr();

    // Numerical pass: accumulate all inputs in their original matrix order
    parallel::for_index(cols, workers, [&](Eigen::Index col, int worker) {
        ColumnBitmap& bitmap = bitmaps[static_cast<std::size_t>(worker)];
        std::size_t   pos    = offsets[static_cast<std::size_t>(col)];

        // Accumulate each stored contribution directly by its global row index
        for (std::size_t input = 0; input < matrices.size(); ++input) {
            const StorageIndex  begin    = outer[input][col];
            const StorageIndex  end      = outer[input][col + 1];
            const StorageIndex* rows_ptr = inner[input];
            const Precision*    vals_ptr = values[input];

            for (StorageIndex i = begin; i < end; ++i) {
                bitmap.add(static_cast<std::size_t>(rows_ptr[i]), vals_ptr[i]);
            }
        }

        // Emit rows in ascending order and reset the hierarchy for reuse
        bitmap.flush([&](std::size_t row) {
            output_rows[pos] = static_cast<StorageIndex>(row);
            output_vals[pos] = bitmap.values[row];
            ++pos;
        });
    });

    return result;
}

} // namespace

/**
 * Reduces one or more square compressed CSC matrices into one global matrix.
 *
 * The public entry point validates the common matrix contract and selects
 * the algorithm. One input is returned by value, two inputs use a direct
 * two-pointer merge, and three or more inputs use hierarchical bitmap
 * accumulation. Inputs remain unchanged, stored zeros are preserved, and
 * duplicate coordinates are summed in input-matrix order.
 *
 * @param matrices Compressed input matrices with identical square dimensions.
 * @return Compressed CSC sum with ascending, unique rows in each column.
 */
SparseMatrix merge_sparse_matrices(const std::vector<SparseMatrix>& matrices) {
    // Validate the common assembly contract before accessing raw CSC data
    logging::error(!matrices.empty(),
        "Sparse merge requires at least one input matrix");

    const Eigen::Index rows = matrices.front().rows();
    const Eigen::Index cols = matrices.front().cols();

    logging::error(rows == cols,
        "Sparse merge requires square matrices");
    for (const SparseMatrix& matrix : matrices) {
        logging::error(matrix.rows() == rows && matrix.cols() == cols,
            "Sparse merge requires identically sized matrices");
        logging::error(matrix.isCompressed(),
            "Sparse merge requires compressed CSC inputs");
    }

    // Select the merge implementation; empty matrices require no accumulation
    if (matrices.size() == 1) {
        return matrices.front();
    }
    if (rows == 0) {
        return SparseMatrix(rows, cols);
    }
    if (matrices.size() == 2) {
        return merge_sparse_pair(matrices[0], matrices[1]);
    }
    return merge_sparse_bitmap(matrices);
}

} // namespace fem::mattools
