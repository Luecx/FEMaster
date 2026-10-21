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
 * compressed partial matrices are combined column by column. Two inputs use
 * a sorted two-pointer merge; three or more use a dense row accumulator with
 * a hierarchical occupancy bitmap. The first pass counts unique structural
 * entries. The second writes the compressed sparse pattern and accumulates
 * numerical contributions in ascending input-worker order.
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
 * @see merge_sparse_matrices
 *
 * @author Finn Eggers
 * @date 28.08.2024
 */

#include "merge_sparse_matrices.h"

#include "../core/logging.h"
#include "../core/parallel.h"
#include "../core/memory.h"
#include "../data/field.h"
#include "../model/element/element.h"

#include <algorithm>
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
 * Estimate the approximate amount of triplets that occur in the model, counting duplicates.
 * This is used to preallocate triplet arrays per thread. Ideally we collect all triplets and then fill it into
 * the matrices
 *
 * @param elements
 * @return
 */
inline uint64_t estimate_triplets(const std::vector<model::ElementPtr>& elements, const SystemDofIds& indices) {
    // total approximate estimated triplets
    uint64_t sum = 0;
    // sparsity of this check, the higher, the more approximate.
    constexpr int sparsity = 16;
    for (size_t i = 0; i < elements.size(); i += sparsity) {
        if (elements[i]) {
            // take the dofs of the first node and assume all other nodes have the same amount of dofs
            // there are no elements with no nodes, so this is safe
            ID id1    = elements[i]->nodes()[0];
            auto dofs = (indices.row(id1).array() >= 0).count();
            sum += elements[i]->n_nodes() * elements[i]->n_nodes() * dofs * dofs * sparsity;
        }
    }

    return sum;
}

/**
 * @brief Flushes accumulated triplets into a sparse matrix.
 *
 * Uses Eigen's setFromTriplets() when the matrix is empty, avoiding the
 * additional sparse-matrix merge performed by insertFromTriplets().
 * Subsequent batches are accumulated using insertFromTriplets().
 *
 * Duplicate entries are summed by Eigen. Empty batches are skipped to
 * avoid unnecessary matrix operations.
 *
 * The triplet buffer is cleared after successful insertion while retaining
 * its allocated capacity for reuse.
 *
 * @param matrix Sparse matrix receiving the accumulated contributions.
 * @param triplets Triplet buffer to insert and clear.
 */
inline void flush_triplets(SparseMatrix& matrix, TripletList& triplets) {
    if (triplets.empty()) {
        return;
    }

    if (matrix.nonZeros() == 0) {
        matrix.setFromTriplets(triplets.begin(), triplets.end());
    } else {
        matrix.insertFromTriplets(triplets.begin(), triplets.end());
    }

    triplets.clear();
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
    // The serial path must also restore linear-algebra threading on return or throw.
    const parallel::ScopedLinearAlgebraThreads threading;

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
    SparseMatrix global_matrix(global_size, global_size);

    // Flush the triplet buffer after approximately sixteen million entries to
    // bound peak temporary memory consumption.
    const size_t MEM_AVAIL  = memory::available_memory_bytes() / 4;
    const size_t BATCH_SIZE = MEM_AVAIL / (sizeof(Triplet));
    const size_t RESERVE    = std::min(BATCH_SIZE, estimate_triplets(elements, indices));

    TripletList triplets{};
    triplets.reserve(RESERVE);

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
        alignas(64) Precision local_matrix_storage[MAX_LOCAL_MATRIX_SIZE * MAX_LOCAL_MATRIX_SIZE];

        // Evaluate the current element matrix into the temporary storage.
        auto local_matrix = compute_local_matrix_with_optional_forces(
            compute_local_matrix,
            element,
            local_matrix_storage,
            nodal_forces
        );

        // Determine the element topology and infer the uniform number of local
        // DOFs associated with each node.
        int num_nodes         = element->n_nodes();
        int local_matrix_size = local_matrix.rows();
        int dofs_per_node     = local_matrix_size / num_nodes;

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
            flush_triplets(global_matrix, triplets);
        }
    }

    // Insert the final incomplete triplet batch after all elements have been
    // processed.
    flush_triplets(global_matrix, triplets);

    return global_matrix;
}


#ifdef _OPENMP

/**
 * @brief Assembles a global sparse matrix using OpenMP and column-wise reduction.
 *
 * Element-local matrices are evaluated in parallel. Every worker owns an
 * independent triplet buffer and partial sparse matrix, so no synchronization
 * is required while element contributions are generated.
 *
 * After compression, merge_sparse_matrices() combines the partial columns
 * using a two-pointer merge or hierarchical bitmap accumulation. It owns
 * structural counting, global CSC allocation and numerical accumulation in
 * input-worker order; this routine retains responsibility for element
 * evaluation, triplet buffering and optional nodal-force reduction.
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
 * @param num_threads Worker count, limited to 128 by the assembly dispatcher.
 * @param min_batch_size Minimum scheduling chunk size for element evaluation.
 *
 * @return Assembled square sparse matrix over all active global DOFs.
 */
template<typename Lambda>
SparseMatrix assemble_matrix_multithreaded(const std::vector<model::ElementPtr>& elements,
                                           const SystemDofIds& indices,
                                           Lambda&& compute_local_matrix,
                                           model::NodeData* nodal_forces,
                                           int num_threads,
                                           std::size_t min_batch_size) {
    // Keep library-level threading suppressed during local assembly and merge,
    // including paths that throw before the final sparse matrix is returned.
    const parallel::ScopedLinearAlgebraThreads threading;

    // Read the requested OpenMP worker count and determine the dimension of the
    // active global system.
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
    // try to use 25% of available memory at MAX across all threads.
    // the initial reserve size is intended to be the batch size at max but fill in approximately twice
    // the estimated amount of triplets. This accounts for different thread speed and inaccuracies of the counting
    const size_t MEM_AVAIL  = memory::available_memory_bytes() / 4;
    const size_t BATCH_SIZE = MEM_AVAIL / (sizeof(Triplet) * num_threads);
    const size_t RESERVE    = std::min(BATCH_SIZE, estimate_triplets(elements, indices) / num_threads * 2);

    // if the system is low on memory, error out. 1024 is very conservative here, this happens when there are
    // less than 16KB per thread available.
    logging::error(BATCH_SIZE > 1024,
        "assemble_matrix_multithreaded: insufficient memory for triplet buffers");

    for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
        thread_triplets[thread_id].reserve(RESERVE);
    }

    // Measure the principal assembly and merge phases independently.
    Timer timer;
    timer.start();

    // Assemble independent elements into their worker-local buffers. The helper
    // captures worker exceptions and rethrows them on the calling thread.
    parallel::for_index(elements.size(), num_threads,
        [&](std::size_t elem_idx, int thread_id) {
            const auto& element = elements[elem_idx];
            if (!element) {
                return;
            }

            // Select the triplet buffer belonging exclusively to the current
            // OpenMP worker.
            auto& local_triplets = thread_triplets[thread_id];

            // Provide aligned fixed-size storage for the local matrix and the
            // element's mapped global DOF identifiers.
            constexpr int MAX_LOCAL_MATRIX_SIZE = 128;
            alignas(64) Precision local_matrix_storage[MAX_LOCAL_MATRIX_SIZE * MAX_LOCAL_MATRIX_SIZE];
            alignas(64) int       global_dofs[MAX_LOCAL_MATRIX_SIZE];

            // Evaluate the current element matrix into the worker-local
            // temporary storage.
            model::NodeData* local_nodal_forces = nodal_forces
                ? &thread_nodal_forces[static_cast<std::size_t>(thread_id)]
                : nullptr;

            auto local_matrix = compute_local_matrix_with_optional_forces(
                compute_local_matrix,
                element,
                local_matrix_storage,
                local_nodal_forces
            );

            // Infer the complete local algebraic dimension and the uniform
            // number of DOFs associated with each node.
            int num_nodes         = element->n_nodes();
            int local_matrix_size = local_matrix.rows();
            int dofs_per_node     = local_matrix_size / num_nodes;
            int local_dof_count   = num_nodes * dofs_per_node;

            // Resolve every local-to-global DOF mapping once per element.
            // Reusing this compact array avoids repeated map lookups during the
            // quadratic local-matrix traversal.
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

                // Skip the complete local row when its associated DOF is
                // inactive.
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

            // Periodically insert large worker-local batches into the
            // associated partial sparse matrix to bound temporary triplet
            // memory.
            if (local_triplets.size() > BATCH_SIZE) {
                flush_triplets(thread_matrices[thread_id], thread_triplets[thread_id]);
            }
        }, min_batch_size);
    // for_index restores the configured counts when it finishes; the following
    // reduction and sparse merge are also explicitly parallel.
    threading.limit();

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
    // matrix. The subsequent reduction reads the compressed CSC arrays
    // directly; its two-input merge also requires sorted row indices.
#pragma omp parallel for num_threads(num_threads) schedule(static, 1)
    for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
        flush_triplets(thread_matrices[thread_id], thread_triplets[thread_id]);
        thread_matrices[thread_id].makeCompressed();
    }

    // Release the outer collection and all reserved triplet capacities after
    // the partial matrices have absorbed their final batches.
    std::vector<TripletList>{}.swap(thread_triplets);

    timer.stop();
    logging::info(true, "Time for last insertion   : ", timer.elapsed(), " ms");
    timer.start();

    // Merge compressed partial matrices while preserving input-worker order
    auto global_matrix = merge_sparse_matrices(thread_matrices);

    timer.stop();
    logging::info(true, "Time for sparse merge     : ", timer.elapsed(), " ms");
    logging::info(true, "Pattern nonzeros          : ", global_matrix.nonZeros());

    return global_matrix;
}
#endif

/**
 * @brief Selects the serial or OpenMP sparse-matrix assembly implementation.
 *
 * When OpenMP support is unavailable at compile time, the serial implementation
 * is used unconditionally. With OpenMP enabled, parallel sparse assembly
 * is selected only when more than one worker thread has been configured.
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
    // Avoid worker-local sparse matrices for small element batches and systems
    constexpr std::size_t min_batch_size = 64;
    const int num_threads = parallel::worker_count(
        elements.size(), std::min(128, static_cast<int>(global_config.max_threads)), min_batch_size
    );

    if (num_threads > 1 && indices.maxCoeff() + 1 > 128) {
        return assemble_matrix_multithreaded(
            elements, indices, compute_local_matrix, nodal_forces, num_threads, min_batch_size
        );
    }

    return assemble_matrix_singlethreaded(elements, indices, compute_local_matrix, nodal_forces);
#endif
}

} } // namespace fem::mattools
