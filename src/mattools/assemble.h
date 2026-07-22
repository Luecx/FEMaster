/**
 * @file assemble.h
 * @brief Implements serial and parallel assembly of sparse finite-element system matrices.
 *
 * The functions in this file assemble a square sparse global matrix from
 * element-local matrix contributions. A caller-provided function computes the
 * local matrix of each element using aligned temporary storage supplied by the
 * assembly routine.
 *
 * The serial implementation accumulates local coefficients in a single triplet
 * buffer and periodically inserts completed batches into the global sparse
 * matrix.
 *
 * When OpenMP is available, the parallel implementation gives each worker an
 * independent triplet buffer and partial sparse matrix. After local assembly,
 * the sorted compressed columns of all partial matrices are combined by a
 * column-wise k-way merge. The first merge pass counts the unique structural
 * entries of every output column. The second pass writes the compressed sparse
 * pattern and simultaneously accumulates all numerical contributions.
 *
 * Each output column is owned by exactly one OpenMP iteration during the final
 * merge. Consequently, the compressed matrix can be filled without locks,
 * atomic operations or intermediate sparse-matrix additions.
 *
 * Eigen and MKL internal threading are disabled while the explicit OpenMP
 * element loop is active to avoid nested parallelism. Their configured thread
 * counts are restored before returning.
 *
 * @see model::Element
 *
 * @author Finn Eggers
 * @date 28.08.2024
 */

#include "../core/logging.h"
#include "../data/field.h"
#include "../model/element/element.h"

#include <algorithm>
#include <limits>
#include <type_traits>

#ifdef _OPENMP
    #include <omp.h>
#endif

namespace fem { namespace mattools {

/**
 * Evaluates an element-local matrix and optionally exposes a nodal force field.
 *
 * Matrix-only assembly callbacks keep the traditional signature
 *
 *     f(element, storage).
 *
 * Nonlinear tangent assembly may instead use
 *
 *     f(element, storage, nodal_forces),
 *
 * where `nodal_forces` is either the caller-owned serial field or the
 * thread-local force field belonging to the active OpenMP worker.
 *
 * @param compute_local_matrix Element-local matrix callback.
 * @param element Current element.
 * @param local_matrix_storage Aligned dense matrix storage.
 * @param nodal_forces Optional nodal force output field.
 * @return Element-local matrix mapped onto the supplied storage.
 */
template<typename Lambda>
auto compute_local_matrix_with_optional_forces(
    Lambda&                 compute_local_matrix,
    const model::ElementPtr& element,
    Precision*              local_matrix_storage,
    model::NodeData*        nodal_forces
) {
    if constexpr (std::is_invocable_v<Lambda&, const model::ElementPtr&, Precision*, model::NodeData&>) {
        logging::error(nodal_forces != nullptr,
            "Matrix assembly callback requires a nodal force output field");

        return compute_local_matrix(element, local_matrix_storage, *nodal_forces);
    } else {
        (void) nodal_forces;
        return compute_local_matrix(element, local_matrix_storage);
    }
}

/**
 * @brief Assembles a global sparse matrix using a single execution thread.
 *
 * The routine evaluates the local matrix of every valid element and maps its
 * local rows and columns to active global system DOFs. Matrix entries involving
 * inactive DOFs are omitted.
 *
 * Local coefficients are collected as sparse triplets. When the triplet buffer
 * exceeds a fixed batch size, the buffered entries are inserted into the global
 * sparse matrix and the buffer is reused. This limits the amount of temporary
 * uncompressed storage required for large models.
 *
 * The local-matrix callback receives the current element and aligned temporary
 * storage owned by this function. It must return an Eigen-compatible matrix
 * view whose dimensions correspond to the element's active local DOFs. If
 * `nodal_forces` is supplied, the callback may additionally receive that field
 * as a third argument and accumulate element force contributions into it.
 *
 * @tparam Lambda Callable type used to compute element-local matrices.
 *
 * @param elements Element collection contributing to the global matrix.
 * @param indices Mapping from node identifiers and local DOFs to active global
 *                system indices. Negative indices identify inactive DOFs.
 * @param compute_local_matrix Callable that computes and returns the local
 *                             matrix of one element.
 * @param nodal_forces Optional nodal force field accumulated by callbacks that
 *                     accept a third argument.
 *
 * @return Assembled square sparse matrix over all active global DOFs.
 */
template<typename Lambda>
SparseMatrix assemble_matrix_singlethreaded(const std::vector<model::ElementPtr>& elements,
                                            const SystemDofIds& indices,
                                            Lambda&& compute_local_matrix,
                                            model::NodeData* nodal_forces = nullptr) {
    // Restrict Eigen to one internal worker because the complete assembly is
    // executed by the calling thread.
    Eigen::setNbThreads(1);

#ifdef USE_MKL
    // Apply the same restriction to MKL so element-level kernels do not create
    // additional worker threads.
    mkl_set_num_threads(1);
#endif

    // The highest active global DOF index determines the dimension of the
    // square assembled matrix. Active indices are assumed to be contiguous and
    // zero-based.
    int global_size = indices.maxCoeff() + 1;

    // Prepare the optional nodal output field before element force accumulation.
    if (nodal_forces) {
        nodal_forces->set_zero();
    }

    // Collect local contributions as triplets before inserting them into the
    // global sparse matrix.
    TripletList triplets;
    SparseMatrix global_matrix(global_size, global_size);

    // Flush the triplet buffer after approximately sixteen million entries to
    // bound peak temporary memory consumption.
    constexpr size_t BATCH_SIZE = 16 * 1024 * 1024;

    // Process all elements in their original storage order.
    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& element = elements[elem_idx];

        // Empty element pointers do not contribute to the assembled matrix.
        if (!element) {
            continue;
        }

        // Provide fixed-size cache-line-aligned storage for the local matrix.
        // The callback must not return a matrix exceeding this capacity.
        constexpr int MAX_LOCAL_MATRIX_SIZE = 128;
        alignas(64) Precision local_matrix_storage[MAX_LOCAL_MATRIX_SIZE * MAX_LOCAL_MATRIX_SIZE]{};

        // Evaluate the current element matrix into the temporary storage.
        auto local_matrix = compute_local_matrix_with_optional_forces(
            compute_local_matrix,
            element,
            local_matrix_storage,
            nodal_forces
        );

        // Determine the element topology and infer the uniform number of local
        // DOFs associated with each node.
        int num_nodes = element->n_nodes();
        int local_matrix_size = local_matrix.rows();
        int dofs_per_node = local_matrix_size / num_nodes;

        // Traverse the local matrix in node-block order and map every
        // coefficient to a global row and column.
        for (int i = 0; i < num_nodes; ++i) {
            for (int j = 0; j < num_nodes; ++j) {
                for (int idof = 0; idof < dofs_per_node; ++idof) {
                    for (int jdof = 0; jdof < dofs_per_node; ++jdof) {
                        // Resolve the active global DOF indices belonging to the
                        // current local row and column.
                        int global_row = indices(element->nodes()[i], idof);
                        int global_col = indices(element->nodes()[j], jdof);

                        // Negative mappings identify DOFs excluded from the
                        // current algebraic system.
                        if(global_row < 0 || global_col < 0) {
                            continue;
                        }

                        // Read the coefficient at the corresponding local
                        // node-and-DOF position.
                        Precision value = local_matrix(i * dofs_per_node + idof, j * dofs_per_node + jdof);

                        // Store the contribution as a sparse triplet. Eigen
                        // combines duplicate global entries during insertion.
                        triplets.emplace_back(global_row, global_col, value);
                    }
                }
            }
        }

        // Periodically insert completed batches into the global matrix to keep
        // the temporary triplet storage bounded.
        if (triplets.size() > BATCH_SIZE) {
            global_matrix.insertFromTriplets(triplets.begin(), triplets.end());

            // Reuse the allocated vector capacity for the next batch.
            triplets.clear();
        }
    }

    // Insert the final incomplete triplet batch after all elements have been
    // processed.
    global_matrix.insertFromTriplets(triplets.begin(), triplets.end());

    return global_matrix;
}


#ifdef _OPENMP

/**
 * @brief Assembles a global sparse matrix using OpenMP and a column-wise k-way merge.
 *
 * Element-local matrices are evaluated in parallel. Every worker owns an
 * independent triplet buffer and partial sparse matrix, so no synchronization
 * is required while element contributions are generated.
 *
 * Once the partial matrices have been compressed, each of their columns
 * contains sorted row indices. The final global matrix is built by merging the
 * corresponding sorted columns of all worker matrices:
 *
 * 1. A first k-way merge counts the number of unique rows in each output
 *    column.
 * 2. A prefix sum converts those counts into compressed-column offsets.
 * 3. A second k-way merge writes each unique row once and accumulates all
 *    numerical contributions belonging to that row.
 *
 * The implementation uses fixed-size arrays for the current iterator position
 * and end position of every worker column. The configured thread count is
 * therefore limited to 128.
 *
 * Every output column is processed by one OpenMP iteration, giving the worker
 * exclusive access to its range in the final compressed sparse arrays. The
 * merge consequently requires neither locks nor atomic additions.
 *
 * @tparam Lambda Callable type used to compute element-local matrices.
 *
 * @param elements Element collection contributing to the global matrix.
 * @param indices Mapping from node identifiers and local DOFs to active global
 *                system indices. Negative indices identify inactive DOFs.
 * @param compute_local_matrix Callable that computes and returns the local
 *                             matrix of one element. When `nodal_forces` is not
 *                             null, it may accept a thread-local nodal force
 *                             field as a third argument.
 * @param nodal_forces Optional global nodal force field. The parallel
 *                     implementation accumulates into thread-local copies and
 *                     reduces them into this field after the element loop.
 *
 * @return Assembled square sparse matrix over all active global DOFs.
 */
template<typename Lambda>
SparseMatrix assemble_matrix_multithreaded(const std::vector<model::ElementPtr>& elements,
                                           const SystemDofIds& indices,
                                           Lambda&& compute_local_matrix,
                                           model::NodeData* nodal_forces = nullptr) {
    // Use the storage-index type selected by Eigen for direct access to the
    // compressed sparse matrix arrays.
    using StorageIndex = typename SparseMatrix::StorageIndex;

    // Disable Eigen's internal parallelism because the outer assembly and final
    // merge are parallelized explicitly with OpenMP.
    Eigen::setNbThreads(1);
#ifdef USE_MKL
    // Prevent MKL routines called during local element evaluation from creating
    // nested worker teams.
    mkl_set_num_threads(1);
#endif

    // Read the requested OpenMP worker count and determine the dimension of the
    // active global system.
    const int num_threads = std::min(128, global_config.max_threads);
    const int global_size = indices.maxCoeff() + 1;

    // Give every worker an independent sparse matrix and triplet buffer. This
    // eliminates shared writes during local element assembly.
    std::vector<SparseMatrix> thread_matrices(num_threads, SparseMatrix(global_size, global_size));
    std::vector<TripletList>  thread_triplets(num_threads);

    // Optional element force output is accumulated into one nodal field per
    // worker. This avoids shared writes at nodes connected to elements processed
    // by different threads.
    std::vector<model::NodeData> thread_nodal_forces;
    if (nodal_forces) {
        nodal_forces->set_zero();
        thread_nodal_forces.reserve(static_cast<std::size_t>(num_threads));

        for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
            thread_nodal_forces.emplace_back(
                nodal_forces->name,
                nodal_forces->domain,
                nodal_forces->rows,
                nodal_forces->components
            );
            thread_nodal_forces.back().set_zero();
        }
    }

    // Reserve a large triplet batch for each worker. Large batches are
    // periodically inserted into the corresponding partial sparse matrix.
    constexpr size_t BATCH_SIZE = 16 * 1024 * 1024;
    for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
        thread_triplets[thread_id].reserve(BATCH_SIZE);
    }

    // Measure the principal assembly and merge phases independently.
    Timer timer;
    timer.start();

    // Distribute element blocks statically among the configured workers. Each
    // iteration writes only into storage owned by its current OpenMP thread.
    #pragma omp parallel for schedule(static, 1024) num_threads(num_threads)
    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& element = elements[elem_idx];

        // Null element pointers do not contribute to the global matrix.
        if (!element) {
            continue;
        }

        // Select the triplet buffer belonging exclusively to the current
        // OpenMP worker.
        const int thread_id = omp_get_thread_num();
        auto& local_triplets = thread_triplets[thread_id];

        // Provide aligned fixed-size storage for the local matrix and the
        // element's mapped global DOF identifiers.
        constexpr int MAX_LOCAL_MATRIX_SIZE = 128;
        alignas(64) Precision local_matrix_storage[MAX_LOCAL_MATRIX_SIZE * MAX_LOCAL_MATRIX_SIZE];
        alignas(64) int       global_dofs[MAX_LOCAL_MATRIX_SIZE];

        // Evaluate the current element matrix into the worker-local temporary
        // storage.
        model::NodeData* local_nodal_forces = nodal_forces
            ? &thread_nodal_forces[static_cast<std::size_t>(thread_id)]
            : nullptr;

        auto local_matrix = compute_local_matrix_with_optional_forces(
            compute_local_matrix,
            element,
            local_matrix_storage,
            local_nodal_forces
        );

        // Infer the complete local algebraic dimension and the uniform number
        // of DOFs associated with each node.
        int num_nodes         = element->n_nodes();
        int local_matrix_size = local_matrix.rows();
        int dofs_per_node     = local_matrix_size / num_nodes;
        int local_dof_count   = num_nodes * dofs_per_node;

        // Resolve every local-to-global DOF mapping once per element. Reusing
        // this compact array avoids repeated map lookups during the quadratic
        // local-matrix traversal.
        for (int node = 0; node < num_nodes; ++node) {
            const ID node_id = element->nodes()[node];

            for (int dof = 0; dof < dofs_per_node; ++dof) {
                global_dofs[node * dofs_per_node + dof] = indices(node_id, dof);
            }
        }

        // Emit all active local matrix entries into the current worker's
        // triplet batch.
        for (int local_row = 0; local_row < local_dof_count; ++local_row) {
            const int global_row = global_dofs[local_row];

            // Skip the complete local row when its associated DOF is inactive.
            if (global_row < 0) {
                continue;
            }

            for (int local_col = 0; local_col < local_dof_count; ++local_col) {
                const int global_col = global_dofs[local_col];

                // Assemble coefficients only when both local indices map to
                // active global DOFs.
                if (global_col >= 0) {
                    local_triplets.emplace_back(global_row, global_col, local_matrix(local_row, local_col));
                }
            }
        }

        // Periodically insert large worker-local batches into the associated
        // partial sparse matrix to bound temporary triplet memory.
        if (local_triplets.size() > BATCH_SIZE) {
            thread_matrices[thread_id].insertFromTriplets(local_triplets.begin(), local_triplets.end());
            local_triplets.clear();
        }
    }

    timer.stop();
    logging::info(true, "Time for parallel loop    : ", timer.elapsed(), " ms");
    timer.start();

    // Reduce the optional thread-local force fields back into the caller-owned
    // field while the partial sparse matrices are still independent.
    if (nodal_forces) {
        const Index value_count = nodal_forces->rows * nodal_forces->components;
        Precision*  target      = nodal_forces->data();

        #pragma omp parallel for num_threads(num_threads) schedule(static, 1024)
        for (Index value = 0; value < value_count; ++value) {
            Precision sum = Precision(0);

            for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
                sum += thread_nodal_forces[static_cast<std::size_t>(thread_id)].data()[value];
            }

            target[value] = sum;
        }

        std::vector<model::NodeData>{}.swap(thread_nodal_forces);

        timer.stop();
        logging::info(true, "Time for force reduction  : ", timer.elapsed(), " ms");
        timer.start();
    }

    // Flush the remaining triplets of every worker and compress each partial
    // matrix. Compression guarantees sorted row indices inside every column,
    // which is required by the subsequent k-way merge.
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
    for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
        thread_matrices[thread_id].insertFromTriplets(thread_triplets[thread_id].begin(),
                                                      thread_triplets[thread_id].end());
        thread_triplets[thread_id].clear();
        thread_matrices[thread_id].makeCompressed();
    }

    // Release the outer collection and all reserved triplet capacities after
    // the partial matrices have absorbed their final batches.
    std::vector<TripletList>{}.swap(thread_triplets);

    timer.stop();
    logging::info(true, "Time for last insertion   : ", timer.elapsed(), " ms");
    timer.start();

    // Store one entry per output column plus the terminal compressed-column
    // offset. During the first pass, column_offsets[col + 1] temporarily holds
    // only the number of unique rows in column col.
    std::vector<Index> column_offsets(global_size + 1, 0);

    // Count unique structural entries by performing an independent k-way merge
    // of every final column.
#pragma omp parallel for num_threads(num_threads) schedule(static, 256)
    for (int col = 0; col < global_size; ++col) {
        // Current compressed-array position for each worker's copy of this
        // column.
        StorageIndex thread_pos[128];

        // One-past-the-end compressed-array position for each worker column.
        StorageIndex thread_end[128];

        // Initialize every worker iterator from the partial matrix's compressed
        // outer-index array.
        for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
            thread_pos[thread_id] = thread_matrices[thread_id].outerIndexPtr()[col];
            thread_end[thread_id] = thread_matrices[thread_id].outerIndexPtr()[col + 1];
        }

        Index unique_count = 0;

        // Repeatedly select the smallest unprocessed row among all worker
        // columns until every iterator reaches its end.
        while (true) {
            StorageIndex next_row = std::numeric_limits<StorageIndex>::max();

            // Inspect the current row of every non-empty worker sequence and
            // determine the globally smallest candidate.
            for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
                if (thread_pos[thread_id] < thread_end[thread_id]) {
                    next_row = std::min(next_row, thread_matrices[thread_id].innerIndexPtr()[thread_pos[thread_id]]);
                }
            }

            // The sentinel remains unchanged only when all worker iterators are
            // exhausted.
            if (next_row == std::numeric_limits<StorageIndex>::max()) {
                break;
            }

            // The selected row contributes exactly one structural nonzero to
            // the final output column, regardless of how many workers contain
            // that row.
            ++unique_count;

            // Advance every worker iterator over all entries matching the
            // selected row. A compressed sparse column normally contains each
            // row once, while the loop also safely handles repeated entries.
            for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
                while (thread_pos[thread_id] < thread_end[thread_id] &&
                       thread_matrices[thread_id].innerIndexPtr()[thread_pos[thread_id]] == next_row) {
                    ++thread_pos[thread_id];
                }
            }
        }

        // Store the structural count temporarily in the following offset slot.
        // Each OpenMP iteration writes to a distinct vector entry.
        column_offsets[col + 1] = unique_count;
    }

    // Convert the per-column structural counts into cumulative compressed-column
    // offsets. After this prefix sum, column_offsets[col] is the beginning of
    // column col in the final inner-index and value arrays.
    for (int col = 0; col < global_size; ++col) {
        column_offsets[col + 1] += column_offsets[col];
    }

    // The terminal offset equals the exact number of nonzero entries required
    // by the final sparse pattern.
    const Index total_nonzeros = column_offsets[global_size];

    // Create the final square matrix and allocate exactly the required amount of
    // compressed sparse storage.
    SparseMatrix global_matrix(global_size, global_size);
    global_matrix.makeCompressed();
    global_matrix.resizeNonZeros(static_cast<Eigen::Index>(total_nonzeros));

    // Obtain direct pointers to Eigen's compressed-column arrays.
    auto* outer_indices = global_matrix.outerIndexPtr();
    auto* inner_indices = global_matrix.innerIndexPtr();
    auto* values        = global_matrix.valuePtr();

    // Repeat the column-wise k-way merge. This pass writes the unique row
    // pattern and simultaneously sums all values belonging to each row.
#pragma omp parallel for num_threads(num_threads) schedule(static, 256)
    for (int col = 0; col < global_size; ++col) {
        // Locate the disjoint output range owned by this column.
        const Index start = column_offsets[col];
        outer_indices[col] = static_cast<StorageIndex>(start);

        StorageIndex thread_pos[128];
        StorageIndex thread_end[128];

        // Initialize the current and terminal positions of every worker column.
        for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
            thread_pos[thread_id] = thread_matrices[thread_id].outerIndexPtr()[col];
            thread_end[thread_id] = thread_matrices[thread_id].outerIndexPtr()[col + 1];
        }

        // Number of unique rows already written into the current final column.
        Index local = 0;

        while (true) {
            StorageIndex next_row = std::numeric_limits<StorageIndex>::max();

            // Find the smallest remaining row across all partial columns.
            for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
                if (thread_pos[thread_id] < thread_end[thread_id]) {
                    next_row = std::min(next_row, thread_matrices[thread_id].innerIndexPtr()[thread_pos[thread_id]]);
                }
            }

            // Stop after all worker sequences have been consumed.
            if (next_row == std::numeric_limits<StorageIndex>::max()) {
                break;
            }

            // Accumulate all numerical contributions associated with the
            // selected row across every worker matrix.
            Precision value = Precision(0);
            for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
                while (thread_pos[thread_id] < thread_end[thread_id] &&
                       thread_matrices[thread_id].innerIndexPtr()[thread_pos[thread_id]] == next_row) {
                    value += thread_matrices[thread_id].valuePtr()[thread_pos[thread_id]];
                    ++thread_pos[thread_id];
                }
            }

            // Write the unique row and its accumulated coefficient into the next
            // entry of this column's preallocated output range.
            const Index offset = start + local++;
            inner_indices[offset] = next_row;
            values[offset]        = value;
        }
    }

    // Complete the compressed-column outer-index array with the terminal offset
    // following the final column.
    outer_indices[global_size] = static_cast<StorageIndex>(total_nonzeros);

    timer.stop();
    logging::info(true, "Time for k-way merge      : ", timer.elapsed(), " ms");
    logging::info(true, "Pattern nonzeros          : ", total_nonzeros);

    // Restore the configured Eigen thread count after explicit OpenMP assembly
    // and merging have completed.
    Eigen::setNbThreads(global_config.max_threads);
#ifdef USE_MKL
    // Restore MKL's configured worker count as well.
    mkl_set_num_threads(global_config.max_threads);
#endif

    return global_matrix;
}
#endif

/**
 * @brief Selects the serial or OpenMP sparse-matrix assembly implementation.
 *
 * When OpenMP support is unavailable at compile time, the serial implementation
 * is used unconditionally. With OpenMP enabled, the parallel k-way assembly is
 * selected only when more than one worker thread has been configured.
 *
 * The element collection, global DOF mapping, local-matrix callback and optional
 * nodal force output are forwarded unchanged to the selected implementation.
 *
 * @tparam Lambda Callable type used to compute element-local matrices.
 *
 * @param elements Element collection contributing to the global matrix.
 * @param indices Mapping from node identifiers and local DOFs to active global
 *                system indices.
 * @param compute_local_matrix Callable that computes and returns one local
 *                             element matrix.
 * @param nodal_forces Optional global nodal force field for callbacks that
 *                     produce internal force contributions.
 *
 * @return Assembled square sparse matrix over all active global DOFs.
 */
template<typename Lambda>
SparseMatrix assemble_matrix(const std::vector<model::ElementPtr>& elements,
                             const SystemDofIds& indices,
                             Lambda&& compute_local_matrix,
                             model::NodeData* nodal_forces = nullptr) {
#ifndef _OPENMP
    // OpenMP support is unavailable, so assembly must remain serial.
    return assemble_matrix_singlethreaded(elements, indices, compute_local_matrix, nodal_forces);
#else
    // Use parallel assembly only when more than one worker is configured.
    if (global_config.max_threads > 1) {
        return assemble_matrix_multithreaded(elements, indices, compute_local_matrix, nodal_forces);
    } else {
        return assemble_matrix_singlethreaded(elements, indices, compute_local_matrix, nodal_forces);
    }
#endif
}

} } // namespace fem::mattools
