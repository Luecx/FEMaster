/**
 * @file assemble.h
 * @brief Assembles sparse finite-element matrices into a fixed shared pattern.
 *
 * The mattools assembly routine obtains a zero-valued compressed sparse matrix
 * from build_sparsity_pattern using element connectivity and the active global
 * DOF map. Pattern construction is independent of local matrix evaluation and
 * is implemented separately in build_sparsity_pattern.cpp.
 *
 * Element callbacks own local matrix evaluation and material-state updates.
 * Assembly evaluates each non-null element once and scatters its coefficients
 * into existing global slots. OpenMP workers share one sparse matrix and use
 * atomic numerical additions, so no triplet buffers, partial matrices or sparse
 * merge are needed. The serial path uses the same algorithm without atomics.
 *
 * Callbacks that also produce nodal forces retain worker-local force fields
 * and a final field reduction. Pattern construction never evaluates callbacks
 * or changes material state. The pattern is rebuilt on every assembly call;
 * caching and reuse remain outside this routine's current responsibility.
 *
 * Eigen and MKL internal threading are suppressed during explicit assembly
 * loops and restored on return or exception by the parallel utilities.
 *
 * @see model::Element
 * @see build_sparsity_pattern
 * @see parallel::for_index
 *
 * @author Finn Eggers
 * @date 28.08.2024
 */

#pragma once

#include "build_sparsity_pattern.h"

#include "../core/logging.h"
#include "../core/parallel.h"
#include "../core/timer.h"
#include "../data/field.h"
#include "../model/element/element.h"

#include <algorithm>
#include <type_traits>

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
 * @brief Assembles a global matrix using a shared connectivity-derived pattern.
 *
 * build_sparsity_pattern supplies the zero-valued compressed matrix without
 * invoking the callback. Its conservative pattern includes all active DOFs of
 * connected nodes, so mixed element families and smaller or empty local matrices
 * can leave explicit structural zeros. Negative DOF identifiers are inactive;
 * active identifiers must be unique. Unused global columns remain empty.
 *
 * Once the compressed structure is allocated, each element callback is evaluated
 * exactly once. Its node-major local matrix contributes A(I_i,I_j) += A_e(i,j)
 * for active global indices I_i and I_j. The local layout uses the first d DOFs
 * per node, with d = local_rows / element_nodes, as in the existing element APIs.
 * Binary search locates each pre-existing coefficient; parallel updates use
 * OpenMP atomics without modifying the structure. With one worker or without
 * OpenMP, the same numerical updates execute without synchronization.
 *
 * The callback receives aligned storage for at most 128 local DOFs and must
 * respect that capacity before returning its matrix view. Returned dimensions
 * are checked before scattering, and sparse-index overflow is checked before
 * coefficient allocation. Null element slots contribute nothing. A mapping
 * without active DOFs produces a zero-dimensional matrix.
 *
 * Optional nodal forces are reset on entry to numerical assembly. Parallel
 * callbacks receive worker-local fields, reduced into the caller's field after
 * assembly; serial callbacks accumulate directly into it. Callback exceptions
 * propagate through parallel::for_index, and internal library thread counts are
 * restored on every exit. Material-state rollback remains the caller's task.
 *
 * @tparam Lambda Element-matrix callback, optionally accepting a force field.
 * @param elements Element connectivity and callback inputs.
 * @param indices Node/DOF to global index map; negative entries are inactive.
 * @param compute_local_matrix Callback returning the local matrix in supplied storage.
 * @param nodal_forces Optional caller-owned nodal force output, reset before accumulation.
 * @return Assembled global matrix with its conservative compressed-column pattern.
 */
template<typename Lambda>
SparseMatrix assemble_matrix(const std::vector<model::ElementPtr>& elements,
                             const SystemDofIds& indices,
                             Lambda&& compute_local_matrix,
                             model::NodeData* nodal_forces = nullptr) {
    const parallel::ScopedLinearAlgebraThreads threading;

    // Bound the element team by useful work. The indexed helper selects one
    // worker without OpenMP or inside an existing OpenMP region. This count
    // sizes the optional force fields and controls numerical accumulation.
    constexpr std::size_t min_batch_size = 64;
    const int num_threads = parallel::worker_count(
        elements.size(), static_cast<int>(global_config.max_threads), min_batch_size);

    // Allocate the zero-valued matrix from connectivity alone. Element
    // evaluation and optional material-state updates begin only afterwards.
    Timer timer;
    timer.start();
    SparseMatrix global_matrix = build_sparsity_pattern(elements, indices);
    threading.limit();
    timer.stop();
    logging::info(true, "Time for sparsity pattern : ", timer.elapsed(), " ms");
    logging::info(true, "Pattern nonzeros          : ", global_matrix.nonZeros());
    timer.start();

    // Matrix coefficients are shared, but callbacks may also scatter nodal
    // forces through the existing non-atomic field API. Keep those fields local
    // to each worker. Serial callbacks write directly to the caller's field.
    std::vector<model::NodeData> thread_nodal_forces;
    if (nodal_forces) {
        nodal_forces->set_zero();
        if (num_threads > 1) {
            thread_nodal_forces.reserve(static_cast<std::size_t>(num_threads));
            for (int worker = 0; worker < num_threads; ++worker) {
                thread_nodal_forces.emplace_back(nodal_forces->name, nodal_forces->domain,
                                                nodal_forces->rows, nodal_forces->components);
                thread_nodal_forces.back().set_zero();
            }
        }
    }

    // The compressed structure stays immutable during all element updates.
    // Only its values are shared between workers; insertion and sparse-matrix
    // reduction are unnecessary because every admissible position exists.
    Precision* values = global_matrix.valuePtr();

    // Evaluate each element once and scatter into existing coefficient slots.
    // Column-major traversal matches Eigen's local and global matrix storage.
    parallel::for_index(elements.size(), num_threads, [&](std::size_t e, int worker) {
        const auto& element = elements[e];
        if (!element) {
            return;
        }
        constexpr int capacity = 128;
        alignas(64) Precision storage[capacity * capacity];
        int global_dofs[capacity];
        model::NodeData* local_nodal_forces = thread_nodal_forces.empty()
            ? nodal_forces : &thread_nodal_forces[worker];
        auto local = compute_local_matrix_with_optional_forces(
            compute_local_matrix, element, storage, local_nodal_forces);
        const int nodes = element->n_nodes();
        const int size  = static_cast<int>(local.rows());
        logging::error(nodes > 0 && size <= capacity && local.cols() == size && size % nodes == 0,
            "Pattern assembly requires a square node-major local matrix of at most 128 DOFs");
        const int dofs_per_node = size / nodes;
        logging::error(dofs_per_node <= indices.cols(),
            "Pattern assembly local DOF layout exceeds the global mapping");
        for (int i = 0; i < size; ++i) {
            global_dofs[i] = indices(element->nodes()[i / dofs_per_node], i % dofs_per_node);
        }
        for (int j = 0; j < size; ++j) {
            const int col = global_dofs[j];
            if (col < 0) {
                continue;
            }
            const auto begin = global_matrix.innerIndexPtr() + global_matrix.outerIndexPtr()[col];
            const auto end   = global_matrix.innerIndexPtr() + global_matrix.outerIndexPtr()[col + 1];
            for (int i = 0; i < size; ++i) {
                const int row = global_dofs[i];
                if (row < 0) {
                    continue;
                }
                // Connectivity guarantees that this row is present in the
                // fixed column. Resolve its coefficient slot without changing
                // Eigen's compressed arrays or invoking sparse insertion.
                const auto      position = std::lower_bound(begin, end, row);
                const Index     offset   = position - global_matrix.innerIndexPtr();
                const Precision value    = local(i, j);
                if (num_threads > 1) {
#ifdef _OPENMP
                    #pragma omp atomic update
#endif
                    values[offset] += value;
                } else {
                    // A single worker needs no synchronization and follows the
                    // same pattern and local-to-global mapping as parallel runs.
                    values[offset] += value;
                }
            }
        }
    }, min_batch_size);
    threading.limit();
    timer.stop();
    logging::info(true, "Time for matrix assembly  : ", timer.elapsed(), " ms");
    timer.start();

    // The matrix is already complete. Only the optional worker-local nodal
    // forces need a reduction; each output value has one writer in this loop.
    if (!thread_nodal_forces.empty()) {
        parallel::for_index(nodal_forces->rows * nodal_forces->components, num_threads, [&](Index i, int) {
            Precision sum = Precision(0);
            for (int worker = 0; worker < num_threads; ++worker) {
                sum += thread_nodal_forces[worker].data()[i];
            }
            nodal_forces->data()[i] = sum;
        });
        timer.stop();
        logging::info(true, "Time for force reduction  : ", timer.elapsed(), " ms");
    }
    return global_matrix;
}

} } // namespace fem::mattools
