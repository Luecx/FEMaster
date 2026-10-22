/**
 * @file build_sparsity_pattern.cpp
 * @brief Builds compressed sparse matrix patterns from finite-element topology.
 *
 * Compact node-to-element incidence supplies each node's neighboring nodes.
 * Their active global DOFs define one sorted row set shared by all columns of
 * that node. Prefix sums allocate the final compressed-column arrays directly,
 * avoiding element matrix evaluation, sparse insertion and structural merging.
 *
 * The returned matrix owns both the pattern and zero-valued coefficient storage.
 * Numerical assembly, material state and nodal forces remain outside this file.
 *
 * @see build_sparsity_pattern
 * @see assemble_matrix
 *
 * @author Finn Eggers
 * @date 22.10.2026
 */

#include "build_sparsity_pattern.h"

#include "../core/logging.h"
#include "../core/parallel.h"
#include "../model/element/element.h"

#include <algorithm>
#include <limits>

namespace fem { namespace mattools {

/**
 * Builds a zero-valued sparse matrix containing every admissible coupling.
 *
 * Two connectivity passes construct compact node-to-element incidence. For
 * each node, the union of connected element nodes is expanded to active global
 * DOFs. Sorting and removing duplicates produces the row set for all columns
 * belonging to that node. A prefix sum allocates the final compressed-column
 * structure, whose numerical values are initialized to zero.
 *
 * All active DOFs of neighboring nodes are included. This is a conservative
 * pattern: mixed element families and callbacks returning smaller or empty
 * matrices may leave explicit structural zeros during subsequent assembly.
 * No element matrix, geometry or material response is evaluated, so stateful
 * callbacks remain entirely under the control of the numerical assembler.
 *
 * Negative DOF indices are inactive. Active indices must be unique; unused
 * global indices produce empty columns. Null elements add no connectivity.
 * An empty map or a map without active DOFs gives a zero-dimensional matrix.
 * Node identifiers and sparse storage limits are validated before indexing or
 * allocation. Worker exceptions propagate through parallel::for_index, and
 * internal Eigen/MKL thread counts are restored on return or exception.
 *
 * @param elements Elements supplying only node counts and connectivity.
 * @param indices Node/DOF to global index map defining the algebraic system.
 * @return Compressed square sparse matrix with zero coefficients and owned storage.
 */
SparseMatrix build_sparsity_pattern(const std::vector<model::ElementPtr>& elements,
                                    const SystemDofIds& indices) {
    using StorageIndex = SparseMatrix::StorageIndex;
    const parallel::ScopedLinearAlgebraThreads threading;

    // Use the same minimum element batch as numerical assembly. Pattern data
    // depend only on connectivity and numbering; worker selection follows the
    // existing global configuration and remains serial without OpenMP.
    constexpr std::size_t min_batch_size = 64;
    const int num_threads = parallel::worker_count(
        elements.size(), static_cast<int>(global_config.max_threads), min_batch_size);
    const int maximum_dof = indices.size() == 0 ? -1 : indices.maxCoeff();
    logging::error(maximum_dof < std::numeric_limits<StorageIndex>::max(),
        "Sparsity pattern dimension exceeds sparse storage index capacity");
    const int        global_size = std::max(0, maximum_dof + 1);
    const std::size_t node_count  = static_cast<std::size_t>(indices.rows());
    SparseMatrix     global_matrix(global_size, global_size);

    // Derive the graph at node level: all active DOFs on one node have the
    // same conservative neighbors. Construct compact node-to-element incidence
    // with counting and prefix sums before determining those neighborhoods.
    // The two connectivity passes use contiguous storage instead of allocating
    // a separate incidence vector for every global algebraic column.
    std::vector<std::size_t> node_offsets(node_count + 1, 0);
    for (const auto& element : elements) {
        if (element) {
            logging::error(element->n_nodes() > 0,
                "Sparsity pattern requires elements with at least one node");
            for (int i = 0; i < element->n_nodes(); ++i) {
                const ID node = element->nodes()[i];
                logging::error(node >= 0 && static_cast<std::size_t>(node) < node_count,
                    "Sparsity pattern element node is outside the global DOF mapping: ", node);
                ++node_offsets[node + 1];
            }
        }
    }
    for (std::size_t node = 0; node < node_count; ++node) {
        node_offsets[node + 1] += node_offsets[node];
    }
    std::vector<std::size_t> incident_elements(node_offsets.back());
    std::vector<std::size_t> next = node_offsets;
    for (std::size_t e = 0; e < elements.size(); ++e) {
        if (elements[e]) {
            for (int i = 0; i < elements[e]->n_nodes(); ++i) {
                incident_elements[next[elements[e]->nodes()[i]]++] = e;
            }
        }
    }
    std::vector<std::size_t>{}.swap(next);

    // Union neighboring node identifiers before expanding them to active global
    // DOFs. The row set is shared by all active columns belonging to this node.
    // Sort the expanded row identifiers as well: global numbering need not
    // follow the order of the connected node identifiers.
    std::vector<std::vector<StorageIndex>> node_rows(node_count);
    parallel::for_index(node_count, num_threads, [&](std::size_t node, int) {
        if (node_offsets[node] == node_offsets[node + 1]) {
            return;
        }
        std::vector<ID> neighbors;
        for (std::size_t p = node_offsets[node]; p < node_offsets[node + 1]; ++p) {
            const auto& element = elements[incident_elements[p]];
            neighbors.insert(neighbors.end(), element->nodes(), element->nodes() + element->n_nodes());
        }
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
        auto& rows = node_rows[node];
        for (ID neighbor : neighbors) {
            for (int dof = 0; dof < indices.cols(); ++dof) {
                const int row = indices(neighbor, dof);
                if (row >= 0) {
                    rows.push_back(row);
                }
            }
        }
        std::sort(rows.begin(), rows.end());
        rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    });
    threading.limit();

    // Resolve each unique active global column to its owning node. The existing
    // system numbering assigns one global index per active nodal DOF; holes may
    // remain empty. Count structural entries in the wider Eigen dimension type
    // before narrowing to the compressed sparse storage index.
    std::vector<ID> column_nodes(global_size, ID(-1));
    for (std::size_t node = 0; node < node_count; ++node) {
        for (int dof = 0; dof < indices.cols(); ++dof) {
            const int col = indices(node, dof);
            if (col >= 0) {
                logging::error(column_nodes[col] < 0,
                    "Sparsity pattern requires unique active global DOF indices");
                column_nodes[col] = static_cast<ID>(node);
            }
        }
    }
    Index nonzeros = 0;
    for (int col = 0; col < global_size; ++col) {
        global_matrix.outerIndexPtr()[col] = static_cast<StorageIndex>(nonzeros);
        if (column_nodes[col] >= 0) {
            nonzeros += static_cast<Index>(node_rows[column_nodes[col]].size());
        }
        logging::error(nonzeros <= std::numeric_limits<StorageIndex>::max(),
            "Sparsity pattern exceeds sparse storage index capacity");
    }
    global_matrix.outerIndexPtr()[global_size] = static_cast<StorageIndex>(nonzeros);
    global_matrix.resizeNonZeros(nonzeros);
    parallel::for_index(global_size, num_threads, [&](int col, int) {
        if (column_nodes[col] < 0) {
            return;
        }
        const auto& rows = node_rows[column_nodes[col]];
        if (rows.empty()) {
            return;
        }
        const auto offset = global_matrix.outerIndexPtr()[col];
        std::copy(rows.begin(), rows.end(), global_matrix.innerIndexPtr() + offset);
        std::fill_n(global_matrix.valuePtr() + offset, rows.size(), Precision(0));
    });
    threading.limit();
    return global_matrix;
}

} } // namespace fem::mattools
