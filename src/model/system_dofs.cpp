/**
 * @file system_dofs.cpp
 * @brief Implements the node-by-component system DOF numbering.
 *
 * The sparsity pattern is built from element connectivity alone. Compact
 * node-to-element incidence supplies each node's neighboring nodes, whose
 * active global DOFs define one sorted row set shared by all columns of that
 * node. Prefix sums allocate the compressed-column structure directly.
 *
 * @see src/model/system_dofs.h
 * @author Finn Eggers
 * @date 19.11.2026
 */

#include "system_dofs.h"

#include "../core/logging.h"
#include "../core/parallel.h"
#include "element/element.h"

#include <algorithm>
#include <limits>

namespace fem::model {

SystemDofs::SystemDofs(Index node_count, Index dof_count)
    : _mask(node_count, dof_count),
      _ids (node_count, dof_count) {
    _mask.fill(false);
    _ids .fill(-1);
}

void SystemDofs::activate(Index node, Index dof) {
    logging::error(!_numerated,
        "SystemDofs: cannot activate DOFs after numbering");
    _mask(node, dof) = true;
}

void SystemDofs::activate(Index node, const Dofs& dofs) {
    logging::error(!_numerated,
        "SystemDofs: cannot activate DOFs after numbering");
    logging::error(width() <= dofs.cols(),
        "SystemDofs: DOF selection is narrower than the numbering");
    for (Index dof = 0; dof < width(); ++dof) {
        _mask(node, dof) = _mask(node, dof) || dofs(dof);
    }
}

void SystemDofs::numerate() {
    logging::error(!_numerated,
        "SystemDofs: numbering already assigned");

    int next = 0;
    for (Index node = 0; node < nodes(); ++node) {
        for (Index dof = 0; dof < width(); ++dof) {
            _ids(node, dof) = _mask(node, dof) ? next++ : -1;
        }
    }

    _size      = next;
    _numerated = true;
}

DynamicVector SystemDofs::reduce(const Field& field) const {
    logging::error(_numerated,
        "SystemDofs: reduce requires a numbering");
    logging::error(field.rows == nodes() && field.components == width(),
        "SystemDofs: field ", field.name, " is ", field.rows, "x", field.components,
        ", numbering is ", nodes(), "x", width());

    DynamicVector vector(_size);
    for (Index node = 0; node < nodes(); ++node) {
        for (Index dof = 0; dof < width(); ++dof) {
            const int id = _ids(node, dof);
            if (id >= 0) vector(id) = field(node, dof);
        }
    }
    return vector;
}

Field SystemDofs::expand(const DynamicVector& vector) const {
    logging::error(_numerated,
        "SystemDofs: expand requires a numbering");
    logging::error(vector.size() == _size,
        "SystemDofs: vector has ", vector.size(), " entries, numbering has ", _size);

    Field field{"EXPANDED_VECTOR", FieldDomain::NODE, nodes(), width()};
    field.set_zero();
    for (Index node = 0; node < nodes(); ++node) {
        for (Index dof = 0; dof < width(); ++dof) {
            const int id = _ids(node, dof);
            if (id >= 0) field(node, dof) = vector(id);
        }
    }
    return field;
}

/**
 * Builds a zero-valued sparse matrix containing every admissible coupling.
 *
 * Every active DOF of a neighboring node is included, yielding a conservative
 * pattern that may contain structural zeros. Null elements are ignored. An
 * empty element collection contributes no entries but preserves the global
 * dimension, which includes the regular-element assembly of POINTMASS-only
 * models whose auxiliary elements are assembled separately.
 *
 * Construction proceeds in five stages:
 *  1. Validate input and count node-to-element incidence.
 *  2. Construct compact incidence arrays using prefix sums.
 *  3. Determine sorted and unique global row indices for each node.
 *  4. Count column entries and allocate the compressed sparse matrix.
 *  5. Populate the CSC row indices and zero-valued coefficients.
 */
SparseMatrix SystemDofs::pattern(const std::vector<ElementPtr>& elements) const {
    using StorageIndex = SparseMatrix::StorageIndex;

    const parallel::ScopedLinearAlgebraThreads threading{};

    constexpr std::size_t min_batch_size = 64;

    logging::error(_numerated,
        "SystemDofs: pattern requires a numbering");

    const std::size_t node_count = static_cast<std::size_t>(nodes());

    // ---------------------------------------------------------------------
    // 1. Validate input and count node-to-element incidence
    // ---------------------------------------------------------------------

    // The numbering, rather than the element collection, determines the
    // matrix dimension
    logging::error(_size < std::numeric_limits<StorageIndex>::max(),
        "Sparsity pattern dimension exceeds sparse storage index capacity");
    const std::size_t global_size = static_cast<std::size_t>(_size);

    // Count how often each node occurs in the element connectivity.
    // Invalid identifiers are recorded without accessing invalid offsets.
    std::vector<std::size_t> node_offsets(node_count + 1, 0);

    bool valid_element_sizes = true;
    bool valid_node_ids      = true;
    ID   invalid_node        = ID(-1);

    for (const auto& element : elements) {
        if (!element) {
            continue;
        }

        const int n_nodes = element->n_nodes();
        if (n_nodes <= 0) {
            valid_element_sizes = false;
            continue;
        }

        const ID* nodes = element->nodes();

        for (int i = 0; i < n_nodes; ++i) {
            const ID node = nodes[i];

            if (node < 0 || static_cast<std::size_t>(node) >= node_count) {
                valid_node_ids = false;
                invalid_node   = node;
                continue;
            }

            ++node_offsets[node + 1];
        }
    }

    // Validate element connectivity before constructing incidence arrays.
    logging::error(valid_element_sizes,
        "Sparsity pattern requires elements with at least one node");
    logging::error(valid_node_ids,
        "Sparsity pattern element node is outside the global DOF mapping: ", invalid_node);

    const int num_threads = parallel::worker_count(
        node_count, static_cast<int>(global_config.max_threads), min_batch_size);

    // ---------------------------------------------------------------------
    // 2. Construct compact node-to-element incidence
    // ---------------------------------------------------------------------

    // Convert incidence counts into offsets using an exclusive prefix sum.
    //
    // All elements incident to node n are stored in the contiguous range:
    //
    //   incident_elements[node_offsets[n] ... node_offsets[n + 1])
    //
    for (std::size_t node = 0; node < node_count; ++node) {
        node_offsets[node + 1] += node_offsets[node];
    }

    std::vector<std::size_t> incident_elements(node_offsets.back());

    // Fill the incidence array using temporary write positions. These
    // positions are discarded immediately after construction.
    {
        std::vector<std::size_t> next = node_offsets;

        for (std::size_t e = 0; e < elements.size(); ++e) {
            const auto& element = elements[e];
            if (!element) {
                continue;
            }

            const int n_nodes = element->n_nodes();
            const ID* nodes   = element->nodes();

            for (int i = 0; i < n_nodes; ++i) {
                incident_elements[next[nodes[i]]++] = e;
            }
        }
    }

    // ---------------------------------------------------------------------
    // 3. Determine sorted global row indices for every mesh node
    // ---------------------------------------------------------------------

    // All active columns belonging to one node share the same conservative
    // row pattern. Neighborhood construction is independent for each node
    // and can therefore be parallelized without synchronization.
    std::vector<std::vector<StorageIndex>> node_rows(node_count);

    parallel::for_index(node_count, num_threads, [&](std::size_t node, int) {
        const std::size_t begin = node_offsets[node];
        const std::size_t end   = node_offsets[node + 1];

        if (begin == end) {
            return;
        }

        // Collect all nodes belonging to incident elements and remove
        // duplicates to obtain the topological neighborhood.
        std::vector<ID> neighbors;

        for (std::size_t p = begin; p < end; ++p) {
            const auto& element = elements[incident_elements[p]];
            const int   n_nodes = element->n_nodes();
            const ID*   nodes   = element->nodes();

            neighbors.insert(neighbors.end(), nodes, nodes + n_nodes);
        }

        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());

        // Expand neighboring nodes into active global DOFs. A second sort
        // is necessary because algebraic numbering need not follow node IDs.
        auto& rows = node_rows[node];

        for (ID neighbor : neighbors) {
            for (int dof = 0; dof < width(); ++dof) {
                const int row = _ids(neighbor, dof);
                if (row >= 0) {
                    rows.push_back(row);
                }
            }
        }

        std::sort(rows.begin(), rows.end());
        rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    });

    // ---------------------------------------------------------------------
    // 4. Count column entries and allocate CSC storage
    // ---------------------------------------------------------------------

    // Reuse the CSC outer-index array as temporary column-length storage.
    // A value of -1 identifies an unused column and also allows duplicate
    // global DOF indices to be detected, including columns with zero rows.
    SparseMatrix global_matrix(global_size, global_size);

    StorageIndex* outer = global_matrix.outerIndexPtr();
    std::fill_n(outer, global_size, StorageIndex(-1));

    int duplicate_column = -1;

    for (std::size_t node = 0; node < node_count; ++node) {
        const auto& rows = node_rows[node];

        for (int dof = 0; dof < width(); ++dof) {
            const int col = _ids(node, dof);
            if (col < 0) {
                continue;
            }

            if (outer[col] >= 0) {
                duplicate_column = col;
            }

            // Each active DOF of this node has the same row count.
            // Row-count capacity is checked before storing the length.
            outer[col] = rows.size() <= static_cast<std::size_t>(
                std::numeric_limits<StorageIndex>::max())
                ? static_cast<StorageIndex>(rows.size())
                : StorageIndex(-2);
        }
    }

    // Determine the total number of structural entries. Negative values
    // identify either unused columns or invalid column lengths.
    Index nonzeros = 0;
    bool valid_column_sizes = true;

    for (int col = 0; col < global_size; ++col) {
        valid_column_sizes &= outer[col] != StorageIndex(-2);

        if (outer[col] >= 0) {
            nonzeros += static_cast<Index>(outer[col]);
        }
    }

    // Validate numbering and storage capacity before converting lengths
    // into CSC offsets or allocating coefficient arrays.
    logging::error(duplicate_column < 0,
        "Sparsity pattern requires unique active global DOF indices: ", duplicate_column);
    logging::error(valid_column_sizes,
        "Sparsity pattern column length exceeds sparse storage index capacity");
    logging::error(nonzeros <= std::numeric_limits<StorageIndex>::max(),
        "Sparsity pattern exceeds sparse storage index capacity");

    // Transform column lengths into an exclusive prefix sum in place.
    // Unused columns contribute zero entries and retain valid empty ranges.
    StorageIndex offset = 0;

    for (int col = 0; col < global_size; ++col) {
        const StorageIndex count = std::max(StorageIndex(0), outer[col]);

        outer[col] = offset;
        offset    += count;
    }

    outer[global_size] = offset;

    // Allocate exactly the required row-index and coefficient storage.
    global_matrix.resizeNonZeros(nonzeros);

    // ---------------------------------------------------------------------
    // 5. Populate CSC row indices and zero-valued coefficients
    // ---------------------------------------------------------------------

    // Each node fills the columns corresponding to its active global DOFs.
    // Unique indices guarantee disjoint storage ranges, so parallel writes
    // require neither locks nor atomic operations.
    parallel::for_index(node_count, num_threads, [&](std::size_t node, int) {
        const auto& rows = node_rows[node];

        if (rows.empty()) {
            return;
        }

        for (int dof = 0; dof < width(); ++dof) {
            const int col = _ids(node, dof);
            if (col < 0) {
                continue;
            }

            const StorageIndex begin = global_matrix.outerIndexPtr()[col];

            std::copy(rows.begin(), rows.end(), global_matrix.innerIndexPtr() + begin);
            std::fill_n(global_matrix.valuePtr() + begin, rows.size(), Precision(0));
        }
    });

    return global_matrix;
}

} // namespace fem::model
