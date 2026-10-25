
/**
 * @file assemble.h
 * @brief Assembles finite-element matrices into a fixed sparse matrix pattern.
 *
 * Global matrix assembly separates sparsity pattern construction from numerical
 * element evaluation. The pattern is either constructed from element connectivity
 * or supplied by the caller for reuse across multiple assembly operations.
 *
 * Each element callback computes a local matrix using aligned temporary storage.
 * Local coefficients are mapped to global DOFs and accumulated directly into
 * the existing compressed-column structure. No sparse insertion, triplet buffers
 * or structural merging are required.
 *
 * Parallel workers share the global coefficient array and use atomic additions
 * for overlapping matrix contributions. Optional nodal forces are accumulated
 * in worker-local fields and reduced after numerical assembly.
 *
 * Pattern ownership, element evaluation and numerical accumulation remain
 * independent. The caller retains responsibility for material-state management
 * and for ensuring that a supplied pattern matches the current system topology.
 *
 * @see build_sparsity_pattern
 * @see model::Element
 * @see parallel::for_index
 *
 * @author Finn Eggers
 * @date 22.10.2026
 */

#pragma once

#include "build_sparsity_pattern.h"

#include "../core/logging.h"
#include "../core/timer.h"
#include "../core/parallel.h"
#include "../core/timer.h"
#include "../data/field.h"
#include "../model/element/element.h"

#include <algorithm>
#include <type_traits>
#include <vector>

namespace fem { namespace mattools {

/**
 * Evaluates an element-local matrix with optional nodal force accumulation.
 *
 * Matrix-only callbacks use the signature
 *
 *     f(element, storage)
 *
 * while callbacks also producing nodal forces use
 *
 *     f(element, storage, nodal_forces).
 *
 * Compile-time dispatch selects the appropriate callback signature without
 * runtime branching. The callback returns a local matrix mapped onto the
 * supplied aligned storage.
 *
 * A force-producing callback requires a valid nodal force field. This
 * requirement is checked by assemble_matrix before numerical assembly.
 *
 * @tparam Lambda Element-local matrix callback type.
 * @param compute_local_matrix Callback evaluating the element matrix.
 * @param element Current finite element.
 * @param local_matrix_storage Aligned storage for the local matrix.
 * @param nodal_forces Optional nodal force accumulation field.
 * @return Element-local matrix mapped onto the supplied storage.
 */
template<typename Lambda>
auto compute_local_matrix_with_optional_forces(
    Lambda&                  compute_local_matrix,
    const model::ElementPtr& element,
    Precision*               local_matrix_storage,
    model::NodeData*         nodal_forces
) {
    if constexpr (std::is_invocable_v<Lambda&, const model::ElementPtr&, Precision*, model::NodeData&>) {
        return compute_local_matrix(element, local_matrix_storage, *nodal_forces);
    } else {
        return compute_local_matrix(element, local_matrix_storage);
    }
}

/**
 * Assembles a global sparse matrix using a fixed compressed-column pattern.
 *
 * The sparsity pattern is either constructed from element connectivity or
 * copied from an existing matrix supplied by the caller. In both cases,
 * numerical coefficients are initialized to zero before assembly.
 *
 * Every non-null element is evaluated exactly once. The local matrix is
 * assembled into the global system according to
 *
 *     K(I_i, I_j) += K_e(i, j),
 *
 * where I_i and I_j are the global indices of the local element DOFs.
 *
 * Local matrices follow node-major ordering. Each element contributes the
 * first d DOFs of every node, with
 *
 *     d = local_matrix.rows() / element.n_nodes().
 *
 * Negative global DOF indices are inactive and contribute nothing.
 *
 * Each global coefficient is located by binary search within its sorted CSC
 * column. Parallel element contributions use atomic additions, while serial
 * execution updates coefficients without synchronization.
 *
 * A supplied pattern must be compressed and compatible with the element
 * connectivity and global DOF mapping. Compatibility is a caller precondition
 * and is not validated in advance. Missing coefficient positions encountered
 * during assembly are nevertheless reported as errors.
 *
 * Optional nodal forces are reset before assembly. Parallel callbacks receive
 * worker-local fields, which are reduced into the caller-owned output field
 * after all elements have been evaluated.
 *
 * Callbacks must not write beyond the supplied capacity of 128 local DOFs.
 * Material-state changes made by callbacks remain the caller's responsibility.
 *
 * @tparam Lambda Element-matrix callback, optionally accepting nodal forces.
 * @param elements Elements supplying connectivity and callback inputs.
 * @param indices Node/DOF to global index map.
 * @param compute_local_matrix Callback evaluating local element matrices.
 * @param nodal_forces Optional caller-owned nodal force output field.
 * @param pattern Optional existing CSC pattern, left unchanged.
 * @return Assembled global matrix with owned compressed sparse storage.
 */
template<typename Lambda>
SparseMatrix assemble_matrix(const std::vector<model::ElementPtr>& elements,
                             const SystemDofIds& indices,
                             Lambda&& compute_local_matrix,
                             model::NodeData* nodal_forces = nullptr,
                             const SparseMatrix* pattern = nullptr) {
    const parallel::ScopedLinearAlgebraThreads threading;

    constexpr std::size_t min_batch_size = 64;
    constexpr int         capacity       = 128;

    // Select the number of workers according to the number of elements.
    // Small models remain serial to avoid unnecessary threading overhead.
    // The same worker count is used for element assembly and force reduction.
    const int num_threads = parallel::worker_count(
        elements.size(), static_cast<int>(global_config.max_threads), min_batch_size);

    // ---------------------------------------------------------------------
    // 1. Validate callback requirements
    // ---------------------------------------------------------------------

    // A callback accepting a nodal force field requires a valid output pointer.
    // This is a compile-time property of the callback, so it only needs to be
    // checked once rather than for every element during numerical assembly.
    constexpr bool requires_forces =
        std::is_invocable_v<Lambda&, const model::ElementPtr&, Precision*, model::NodeData&>;

    logging::error(!requires_forces || nodal_forces != nullptr,
        "Matrix assembly callback requires a nodal force output field");
    logging::error(pattern == nullptr || pattern->isCompressed(),
        "Matrix assembly callback requires a compressed sparsity pattern");

    // ---------------------------------------------------------------------
    // 2. Prepare the global sparse matrix
    // ---------------------------------------------------------------------

    Timer timer;
    timer.start();

    // Construct the sparsity pattern from connectivity unless the caller
    // already provides one. An existing pattern is copied because the
    // returned matrix must own its coefficients independently.
    //
    // Reusing a pattern avoids the topological neighborhood construction,
    // sorting and structural deduplication performed by build_sparsity_pattern.
    //
    // No connectivity validation is repeated here. A supplied pattern is
    // assumed to match the current element topology and global DOF numbering.
    SparseMatrix global_matrix = pattern
        ? *pattern
        : build_sparsity_pattern(elements, indices);

    // A copied matrix may contain coefficients from an earlier assembly.
    // Reset its numerical values while retaining the complete CSC structure.
    // Newly constructed patterns already contain zero-valued coefficients.
    if (pattern) {
        std::fill_n(global_matrix.valuePtr(), global_matrix.nonZeros(), Precision(0));
    }

    threading.limit();
    timer.stop();

    logging::info(true, "Time for pattern setup     : ", timer.elapsed(), " ms");
    logging::info(true, "Pattern nonzeros           : ", global_matrix.nonZeros());

    // ---------------------------------------------------------------------
    // 3. Prepare optional nodal force accumulation
    // ---------------------------------------------------------------------

    // Matrix coefficients support atomic updates, but the existing nodal
    // force field interface does not. Parallel callbacks therefore cannot
    // safely accumulate into the same force field.
    //
    // With one worker, the callback writes directly to nodal_forces.
    // With multiple workers, each worker receives a private field having
    // the same dimensions and domain as the caller-owned output.
    std::vector<model::NodeData> thread_nodal_forces;

    if (nodal_forces) {
        nodal_forces->set_zero();

        if (num_threads > 1) {
            thread_nodal_forces.reserve(num_threads);

            for (int worker = 0; worker < num_threads; ++worker) {
                thread_nodal_forces.emplace_back(
                    nodal_forces->name,
                    nodal_forces->domain,
                    nodal_forces->rows,
                    nodal_forces->components);

                thread_nodal_forces.back().set_zero();
            }
        }
    }

    // Obtain direct access to Eigen's compressed-column storage:
    //
    //   outer[j] ... outer[j + 1]   : stored entries of column j
    //   inner[k]                    : global row index of entry k
    //   values[k]                   : numerical coefficient of entry k
    //
    // The column pointers and row indices remain immutable throughout
    // numerical assembly. Only the numerical coefficient array is modified.
    //
    // Since each column stores sorted row indices, binary search can locate
    // any admissible matrix coefficient without inserting new entries.
    const auto* outer = global_matrix.outerIndexPtr();
    const auto* inner = global_matrix.innerIndexPtr();
    Precision*  values = global_matrix.valuePtr();

    // ---------------------------------------------------------------------
    // 4. Evaluate element matrices and assemble global coefficients
    // ---------------------------------------------------------------------

    timer.start();

    // Elements are evaluated independently. Each worker computes one local
    // element matrix at a time and scatters its coefficients into the shared
    // global CSC storage.
    //
    // The sparsity pattern remains fixed. Different workers may update the
    // same coefficient, but never modify the structure of the sparse matrix.
    parallel::for_index(elements.size(), num_threads, [&](std::size_t e, int worker) {
        const auto& element = elements[e];

        if (!element) {
            return;
        }

        // Provide aligned temporary storage for the element matrix.
        // The callback may return a matrix view of any supported dimension
        // up to capacity, but must never write beyond this storage.
        //
        // global_dofs maps each local matrix row or column to its global
        // algebraic index. Negative indices identify inactive DOFs.
        alignas(64) Precision storage[capacity * capacity];
        int global_dofs[capacity];

        // Select the appropriate nodal force destination.
        //
        // Serial assembly writes directly into the caller's output field.
        // Parallel assembly uses the field belonging to the current worker,
        // avoiding races between element callbacks that share mesh nodes.
        model::NodeData* local_nodal_forces = thread_nodal_forces.empty()
            ? nodal_forces
            : &thread_nodal_forces[worker];

        // Evaluate the local matrix exactly once. The callback may also
        // update material state and accumulate local nodal forces.
        //
        // Neither the sparsity pattern nor the numerical assembly algorithm
        // depends on the element formulation used by this callback.
        auto local = compute_local_matrix_with_optional_forces(
            compute_local_matrix, element, storage, local_nodal_forces);

        const int nodes = element->n_nodes();
        const int size  = static_cast<int>(local.rows());

        // Validate the node-major local matrix layout before constructing
        // the local-to-global index mapping.
        //
        // The matrix must be square, fit the supplied storage and contain
        // an equal number of local DOFs for every element node.
        // A zero-dimensional matrix is valid and contributes no coefficients.
        logging::error(nodes > 0 && size >= 0 && size <= capacity
                    && local.cols() == size && size % nodes == 0,
            "Matrix assembly requires a square node-major local matrix of at most 128 DOFs");

        const int dofs_per_node = size / nodes;

        logging::error(dofs_per_node <= indices.cols(),
            "Matrix assembly local DOF layout exceeds the global mapping");

        // Construct the local-to-global mapping in node-major order:
        //
        //   local node 0: DOF 0, DOF 1, ..., DOF d-1
        //   local node 1: DOF 0, DOF 1, ..., DOF d-1
        //   ...
        //
        // Grouping the mapping by node avoids repeated division and modulo
        // operations when retrieving the global index of each local DOF.
        int k = 0;

        for (int node = 0; node < nodes; ++node) {
            const ID node_id = element->nodes()[node];

            for (int dof = 0; dof < dofs_per_node; ++dof) {
                global_dofs[k++] = indices(node_id, dof);
            }
        }

        // Traverse the local matrix column by column. This follows Eigen's
        // column-major storage convention and reuses the same global CSC
        // column range for all entries belonging to one local column.
        //
        // Only active global DOFs contribute to the assembled system.
        for (int j = 0; j < size; ++j) {
            const int col = global_dofs[j];

            if (col < 0) {
                continue;
            }

            // Resolve the contiguous CSC range belonging to this column.
            // Row indices in this range are sorted in ascending order.
            const auto* begin = inner + outer[col];
            const auto* end   = inner + outer[col + 1];

            for (int i = 0; i < size; ++i) {
                const int row = global_dofs[i];

                if (row < 0) {
                    continue;
                }

                // Find the coefficient K(row, col) within the fixed column.
                //
                // The conservative connectivity-derived pattern guarantees
                // that every admissible element coupling already exists.
                // A supplied pattern must satisfy the same requirement.
                //
                // Retain this runtime check even when an existing pattern is
                // trusted. It detects missing entries before an invalid
                // coefficient offset can be used for numerical accumulation.
                const auto* position = std::lower_bound(begin, end, row);

                logging::error(position != end && *position == row,
                    "Matrix assembly pattern is missing coupling (", row, ", ", col, ")");

                // The position in the global row-index array is also the
                // position of the corresponding numerical coefficient.
                // Both arrays have identical ordering and entry counts.
                const Index     offset = position - inner;
                const Precision value  = local(i, j);

                // Multiple elements may contribute to the same global
                // coefficient. Parallel assembly therefore requires atomic
                // addition to prevent lost updates.
                //
                // With a single worker, synchronization is unnecessary.
                // Both execution paths use the same local-to-global mapping
                // and accumulate into the same fixed CSC structure.
                if (num_threads > 1) {
#ifdef _OPENMP
                    #pragma omp atomic update
#endif
                    values[offset] += value;
                } else {
                    values[offset] += value;
                }
            }
        }
    }, min_batch_size);

    threading.limit();
    timer.stop();

    logging::info(true, "Time for matrix assembly   : ", timer.elapsed(), " ms");

    // ---------------------------------------------------------------------
    // 5. Reduce optional worker-local nodal forces
    // ---------------------------------------------------------------------

    // Numerical matrix assembly is complete at this point.
    //
    // If worker-local force fields were created, their contributions must
    // now be summed into the caller-owned nodal force field:
    //
    //   f_global[i] = sum_worker f_worker[i]
    //
    // The reduction is parallelized over output entries rather than workers.
    // Each output entry is therefore written by exactly one worker and does
    // not require atomic operations or additional synchronization.
    //
    // Serial assembly skips this phase because callbacks already accumulated
    // directly into nodal_forces.
    if (!thread_nodal_forces.empty()) {
        timer.start();

        const Index count = static_cast<Index>(nodal_forces->rows) * nodal_forces->components;
        Precision*  output = nodal_forces->data();

        parallel::for_index(count, num_threads, [&](Index i, int) {
            Precision sum = Precision(0);

            for (int worker = 0; worker < num_threads; ++worker) {
                sum += thread_nodal_forces[worker].data()[i];
            }

            output[i] = sum;
        });

        timer.stop();

        logging::info(true, "Time for force reduction   : ", timer.elapsed(), " ms");
    }

    return global_matrix;
}

} } // namespace fem::mattools
