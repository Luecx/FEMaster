/******************************************************************************
 * @file assemble.cpp
 * @brief Provides the implementation of the assemble_matrix function for system
 * matrices in FEM, using user-defined lambda functions.
 *
 * The assembly process handles the global sparse matrix construction by iterating
 * over a vector of element pointers and calling the lambda for local matrix
 * contributions, which are then assembled into the global matrix.
 *
 * This implementation uses multi-threading where available and also estimates
 * the number of non-zero elements (nnz) to improve memory allocation for the
 * sparse matrix.
 *
 * @author Created by Finn Eggers
 * @date Created on 28.08.2024
 ******************************************************************************/

namespace fem { namespace mattools {

/**
 * @brief Assembles the global sparse matrix from local element matrices using a single thread.
 * This function uses a single vector for the triplets and a single matrix.
 *
 * @tparam Lambda A callable that computes the local matrix for an element.
 * @param elements A vector of pointers to the finite elements in the model.
 * @param indices A mapping from node and DOF to global indices.
 * @param compute_local_matrix A lambda function to compute the local element matrix.
 * @return The assembled global sparse matrix.
 */
template<typename Lambda>
SparseMatrix assemble_matrix_singlethreaded(const std::vector<model::ElementPtr>& elements,
                                            const SystemDofIds& indices,
                                            Lambda&& compute_local_matrix) {
    // Ensure Eigen and MKL use only one thread
    Eigen::setNbThreads(1);
    mkl_set_num_threads(1);

    // Determine the size of the global matrix
    int global_size = indices.maxCoeff() + 1;

    // Prepare a single triplet list and a sparse matrix
    std::vector<Eigen::Triplet<double>> triplets;
    SparseMatrix global_matrix(global_size, global_size);

    // Define batch size for squashing the buffer
    constexpr size_t BATCH_SIZE = 16 * 1024 * 1024;

    // Iterate over all elements
    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& element = elements[elem_idx];
        if (!element) {
            continue;  // Skip null elements
        }

        // Local storage for the element's local matrix
        constexpr int MAX_LOCAL_MATRIX_SIZE = 128;  // Maximum size of local matrices
        alignas(64) double local_matrix_storage[MAX_LOCAL_MATRIX_SIZE * MAX_LOCAL_MATRIX_SIZE]{};

        // Compute the local matrix for the current element
        auto local_matrix = compute_local_matrix(element, local_matrix_storage);

        // Get element-specific data: number of nodes and DOFs per node
        int num_nodes = element->n_nodes();
        int local_matrix_size = local_matrix.rows();
        int dofs_per_node = local_matrix_size / num_nodes;

        // Assemble the local matrix into the triplet list
        for (int i = 0; i < num_nodes; ++i) {
            for (int j = 0; j < num_nodes; ++j) {
                for (int idof = 0; idof < dofs_per_node; ++idof) {
                    for (int jdof = 0; jdof < dofs_per_node; ++jdof) {
                        int global_row = indices(element->nodes()[i], idof);
                        int global_col = indices(element->nodes()[j], jdof);
                        double value = local_matrix(i * dofs_per_node + idof, j * dofs_per_node + jdof);

                        // Add the computed value to the triplet list
                        triplets.emplace_back(global_row, global_col, value);
                    }
                }
            }
        }

        // If triplet buffer exceeds batch size, insert into the global matrix and clear the buffer
        if (triplets.size() > BATCH_SIZE) {
            global_matrix.insertFromTriplets(triplets.begin(), triplets.end());
            triplets.clear();  // Clear the triplet buffer after squashing
        }
    }

    // Insert any remaining triplets after the loop
    global_matrix.insertFromTriplets(triplets.begin(), triplets.end());

    return global_matrix;
}


/**
 * @brief Assembles the global sparse matrix from local element matrices using multiple threads.
 * Each thread computes local matrices and assembles them into its own global matrix.
 * At the end, the matrices from all threads are added together to form the final global matrix.
 *
 * @tparam Lambda A callable that computes the local matrix for an element.
 * @param elements A vector of pointers to the finite elements in the model.
 * @param indices A mapping from node and DOF to global indices.
 * @param compute_local_matrix A lambda function to compute the local element matrix.
 * @return The assembled global sparse matrix.
 */
template<typename Lambda>
SparseMatrix assemble_matrix_multithreaded(const std::vector<model::ElementPtr>& elements,
                                           const SystemDofIds& indices,
                                           Lambda&& compute_local_matrix) {
    // Ensure Eigen and MKL use only one thread each
    Eigen::setNbThreads(1);
    mkl_set_num_threads(1);

    // Determine the number of threads for parallel execution
    int num_threads = global_config.max_threads;  // Assuming global_config.max_threads is defined

    // Determine the size of the global matrix
    int global_size = indices.maxCoeff() + 1;

    // Prepare a vector of per-thread sparse matrices
    std::vector<SparseMatrix> thread_matrices(num_threads, SparseMatrix(global_size, global_size));

    // Each thread will collect triplets for its own matrix
    std::vector<std::vector<Eigen::Triplet<double>>> thread_triplets(num_threads);

    // Define batch size for squashing the buffer
    constexpr size_t BATCH_SIZE = 16 * 1024 * 1024;

    // Main parallel loop over elements
    #pragma omp parallel for schedule(static, 1024) num_threads(num_threads)
    for (size_t elem_idx = 0; elem_idx < elements.size(); ++elem_idx) {
        const auto& element = elements[elem_idx];
        if (!element) {
            continue;  // Skip null elements
        }

        // Get the current thread's ID
        int thread_id = omp_get_thread_num();

        // Get the triplets vector for this thread
        auto& local_triplets = thread_triplets[thread_id];

        // Local storage for the element's local matrix
        constexpr int MAX_LOCAL_MATRIX_SIZE = 128;  // Maximum size of local matrices
        alignas(64) double local_matrix_storage[MAX_LOCAL_MATRIX_SIZE * MAX_LOCAL_MATRIX_SIZE]{};

        // Compute the local matrix for the current element
        auto local_matrix = compute_local_matrix(element, local_matrix_storage);

        // Get element-specific data: number of nodes and DOFs per node
        int num_nodes = element->n_nodes();
        int local_matrix_size = local_matrix.rows();
        int dofs_per_node = local_matrix_size / num_nodes;

        // Assemble the local matrix into the thread's triplet list
        for (int i = 0; i < num_nodes; ++i) {
            for (int j = 0; j < num_nodes; ++j) {
                for (int idof = 0; idof < dofs_per_node; ++idof) {
                    for (int jdof = 0; jdof < dofs_per_node; ++jdof) {
                        int global_row = indices(element->nodes()[i], idof);
                        int global_col = indices(element->nodes()[j], jdof);
                        double value = local_matrix(i * dofs_per_node + idof, j * dofs_per_node + jdof);

                        // Add the computed value to the thread-local triplet list
                        local_triplets.emplace_back(global_row, global_col, value);
                    }
                }
            }
        }

        // If triplet buffer exceeds batch size, insert into the thread-local matrix and clear the buffer
        if (local_triplets.size() > BATCH_SIZE) {
            thread_matrices[thread_id].insertFromTriplets(local_triplets.begin(), local_triplets.end());
            local_triplets.clear();  // Clear the triplet buffer after squashing
        }
    }

#pragma omp parallel for num_threads(num_threads)
    for (int thread_id = 0; thread_id < num_threads; ++thread_id) {
        thread_matrices[thread_id].insertFromTriplets(thread_triplets[thread_id].begin(),
                                                      thread_triplets[thread_id].end());
    }

    for(int i = 1; i < num_threads; i++) {
        thread_matrices[0] += thread_matrices[i];
    }
    // The final result is in thread_matrices[0]
    return thread_matrices[0];
}


template<typename Lambda>
SparseMatrix assemble_matrix(const std::vector<model::ElementPtr>& elements,
                             const SystemDofIds& indices,
                             Lambda&& compute_local_matrix) {
    if (global_config.max_threads > 1) {
        return assemble_matrix_multithreaded(elements, indices, compute_local_matrix);
    } else {
        return assemble_matrix_singlethreaded(elements, indices, compute_local_matrix);
    }
}

} } // namespace fem::mattools
